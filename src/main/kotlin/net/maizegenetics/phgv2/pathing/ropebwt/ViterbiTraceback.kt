package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.ceil
import kotlin.math.sqrt

/**
 * Runs a Viterbi recursion and traces its best path back, holding backpointers for a bounded number of
 * positions at a time.
 *
 * A Viterbi traceback needs, for every position after the first and every state, the predecessor that
 * state's best path came from: `positionCount * stateCount` entries. With 82 founders a diploid path
 * has 6,724 states, so that passes 2^31 -- the largest array the JVM allows -- at about 319,000
 * positions, and is 8.6 GB well before then.
 *
 * When the table fits in [maxStoredBackPointers] entries it is kept whole and nothing changes. When it
 * does not, the positions are split into blocks. The first pass keeps only the score vector at the
 * start of each block (a *checkpoint*) and discards the backpointers. The traceback then works from
 * the last block to the first: it restores the block's checkpoint, reruns [step] over the block to
 * rebuild that block's backpointers, and follows them back to the block's start, which gives the state
 * to begin the previous block from.
 *
 * The rerun repeats exactly the arithmetic of the first pass, from exactly the same starting scores, so
 * the path -- including how ties fall -- is identical to keeping the whole table. The cost is one more
 * call of [step] per position, about twice the time. Memory is one block of backpointers plus one
 * checkpoint per block. A block is as long as the budget allows, and at least `sqrt(positionCount)`
 * positions so that the checkpoints stay small too: at a million positions and 6,724 states the
 * minimum is about 27 MB of backpointers and 54 MB of checkpoints, against 27 GB for the whole table.
 *
 * [step] may be called more than once for a position, so it and anything it reads -- the emission in
 * particular -- must depend only on the position and the scores it is given.
 *
 * @param stateCount the number of states
 * @param positionCount the number of positions
 * @param maxStoredBackPointers the most backpointer entries to hold at once; see [ViterbiHMM]
 * @param initialize fills the scores at position 0
 * @param step fills `current` with the best score into each state at `position` from the `previous`
 *   scores, and `backPointers` with the predecessor state each one came from
 * @return the best state at each position, and that path's score
 */
internal fun traceBestPath(
    stateCount: Int,
    positionCount: Int,
    maxStoredBackPointers: Long,
    initialize: (scores: DoubleArray) -> Unit,
    step: (position: Int, previous: DoubleArray, current: DoubleArray, backPointers: IntArray) -> Unit
): Pair<IntArray, Double> {
    require(stateCount > 0) { "State count must be positive" }
    require(positionCount > 0) { "Position count must be positive" }

    var previous = DoubleArray(stateCount)
    var current = DoubleArray(stateCount)
    initialize(previous)

    // Each block covers the steps into positions (start, start + blockLength], so the steps into
    // positions 1 until positionCount fall into stepCount / blockLength blocks, rounded up.
    val stepCount = positionCount - 1
    val blockLength = backPointerBlockLength(stateCount, stepCount, maxStoredBackPointers)
    val blockCount = if (stepCount == 0) 0 else (stepCount + blockLength - 1) / blockLength
    val singleBlock = blockCount <= 1
    // rows[i] holds the backpointers for the step into position (block start + 1 + i)
    val rows = Array(minOf(blockLength, stepCount)) { IntArray(stateCount) }
    val checkpoints = arrayOfNulls<DoubleArray>(if (singleBlock) 0 else blockCount)

    for (position in 1 until positionCount) {
        val offset = position - 1
        if (!singleBlock && offset % blockLength == 0) checkpoints[offset / blockLength] = previous.copyOf()
        // With one block every row is kept; otherwise row 0 is scratch, since this pass only needs
        // the scores.
        step(position, previous, current, if (singleBlock) rows[offset] else rows[0])
        val swap = previous; previous = current; current = swap
    }

    var bestFinalState = 0
    var bestFinalScore = previous[0]
    for (state in 1 until stateCount) {
        if (previous[state] > bestFinalScore) {
            bestFinalScore = previous[state]; bestFinalState = state
        }
    }

    val best = IntArray(positionCount)
    best[positionCount - 1] = bestFinalState
    if (singleBlock) {
        for (position in positionCount - 1 downTo 1) best[position - 1] = rows[position - 1][best[position]]
        return best to bestFinalScore
    }

    for (block in blockCount - 1 downTo 0) {
        val start = block * blockLength
        val end = minOf(start + blockLength, positionCount - 1)
        var blockPrevious = checkpoints[block]!!
        checkpoints[block] = null
        var blockCurrent = current
        for (position in start + 1..end) {
            step(position, blockPrevious, blockCurrent, rows[position - start - 1])
            val swap = blockPrevious; blockPrevious = blockCurrent; blockCurrent = swap
        }
        current = blockCurrent
        for (position in end downTo start + 1) best[position - 1] = rows[position - start - 1][best[position]]
    }
    return best to bestFinalScore
}

/**
 * How many steps' backpointers [traceBestPath] holds at once: all of them when they fit in
 * [maxStoredBackPointers], otherwise as many as fit, but never fewer than `sqrt(stepCount)`.
 */
internal fun backPointerBlockLength(stateCount: Int, stepCount: Int, maxStoredBackPointers: Long): Int {
    if (stepCount <= 1) return 1
    if (stepCount.toLong() * stateCount <= maxStoredBackPointers) return stepCount
    val withinBudget = (maxStoredBackPointers / stateCount).coerceAtMost(stepCount.toLong()).toInt()
    val minimum = ceil(sqrt(stepCount.toDouble())).toInt()
    return maxOf(withinBudget, minimum, 1)
}

/**
 * A default for [ViterbiHMM.maxStoredBackPointers]: a quarter of the maximum heap, shared among
 * [concurrentPaths] paths found at the same time, at four bytes an entry.
 */
fun defaultBackPointerBudget(concurrentPaths: Int = 1): Long =
    Runtime.getRuntime().maxMemory() / 4 / 4 / concurrentPaths.coerceAtLeast(1)
