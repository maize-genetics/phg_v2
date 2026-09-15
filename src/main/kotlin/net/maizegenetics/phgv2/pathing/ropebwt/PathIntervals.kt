package net.maizegenetics.phgv2.pathing.ropebwt

import net.maizegenetics.phgv2.utils.Position

/**
 * A run of consecutive bins that share the same call, expressed in 0-based half-open
 * [start, end) reference coordinates per the BED specification.
 *
 * @param T the type of the call. Haploid paths use the parent name; diploid paths use an
 *   ordered (parent1, parent2) pair.
 */
data class PathInterval<T>(val contig: String, val start: Int, val end: Int, val call: T)

/**
 * Converts an ordered, per-bin path into reference-coordinate intervals.
 *
 * The cut between bin i and bin i+1 is placed at the midpoint of the two bins,
 * `((pos[i] + pos[i+1]) / 2) * binSize`. Each cut is the exclusive end of the interval before it
 * and the inclusive start of the interval after it, so the returned intervals are gapless and
 * non-overlapping. The first interval starts at 0 and the last ends at `pos[last] * binSize`.
 *
 * Intervals covering no reference bases are not emitted. Only the first run can be empty: it is
 * forced to start at 0 rather than at a cut, so when bins 0 and 1 are both present with different
 * calls the cut between them lands at `(0 + 1) / 2 * binSize` = 0 and the run spans `[0, 0)`.
 * Coordinates stay 0-based half-open throughout; a zero-length record is simply malformed in that
 * convention, and its bases are covered by the run that follows. Downstream consumers reject it:
 * `BedToVcf` converts a record to 1-based inclusive as `(start + 1, end)`, and an `end` of 0 has
 * no valid 1-based inclusive form, so it throws. A `start` of 0 is fine and always was.
 *
 * Note this drops bin 0's call rather than giving it any sequence. That follows from the midpoint
 * rule under integer division and is worth revisiting: an alternative is to floor the first cut at
 * `binSize` so bin 0 keeps `[0, binSize)`. Both are defensible; this one changes no existing
 * interval boundary.
 *
 * When [mergeAdjacent] is true, consecutive bins whose calls are equal are emitted as a single
 * interval and the cuts between them are dropped. Equality is whatever `==` means for [T], so for
 * diploid paths (lineA:0, lineB:0) and (lineB:0, lineA:0) are different calls and are not merged,
 * preserving phase.
 *
 * @param path (position, call) pairs for a single contig, ordered by bin position.
 * @param binSize the bin size used to create the ps4g file.
 * @param mergeAdjacent if false, every bin is emitted as its own interval.
 * @return the intervals in path order, or an empty list if [path] is empty.
 */
fun <T> pathToIntervals(
    path: List<Pair<Position, T>>,
    binSize: Int,
    mergeAdjacent: Boolean = true
): List<PathInterval<T>> {

    if (path.isEmpty()) return emptyList()

    //the cut between bin index and bin index + 1: exclusive end of the former, inclusive start of the latter
    fun cut(index: Int) = (path[index].first.position + path[index + 1].first.position) / 2 * binSize

    val intervals = mutableListOf<PathInterval<T>>()
    var runStart = 0

    for (index in path.indices) {
        val isLast = index == path.lastIndex
        val endsRun = isLast || !mergeAdjacent || path[index].second != path[index + 1].second
        if (!endsRun) continue

        val start = if (runStart == 0) 0 else cut(runStart - 1)
        val end = if (isLast) path[index].first.position * binSize else cut(index)
        if (end > start) {
            intervals.add(PathInterval(path[runStart].first.contig, start, end, path[runStart].second))
        }
        runStart = index + 1
    }

    return intervals
}
