package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.exp
import kotlin.math.expm1
import kotlin.math.ln
import kotlin.math.ln1p

/**
 * Per-step transition log probabilities for observations that are not evenly spaced.
 *
 * A fixed per-observation switch probability makes the penalty for changing founders depend on how
 * many markers there are rather than on how much sequence lies between them: two markers 200 bp apart
 * and two 2 Mb apart would cost the same. Where the true positions are known, as they are when the
 * observations are VCF sites, the penalty can be charged per base instead.
 *
 * Switching is modelled as a Poisson process, so over a distance `d`
 *
 *     P(stay on the same founder) = exp(-r * d),   r = -ln(1 - probSwitch) / referenceDistance
 *
 * with [probSwitch] the probability of at least one switch across [referenceDistance].
 *
 * ## Units, and a genetic map later
 *
 * [distances] carries whatever unit [referenceDistance] is quoted in. Physical base pairs are the
 * obvious choice, but centimorgans from a genetic map would substitute directly -- recombination rate
 * varies a great deal along a chromosome, and a map captures that where a single rate cannot. Nothing
 * below or downstream would change; only what the caller measures.
 *
 * One caveat if that is ever done. `impute-path-from-ps4g`'s default, `--prob-switch` 1e-16 per Mb,
 * chosen on the maize simulated corpus, is some thirteen orders of magnitude less recombination than
 * real maize (roughly 0.005 per Mb). That gap is deliberate smoothing against noisy emissions -- each
 * read in a bin counts as independent evidence, which overstates it at high coverage -- not an error,
 * so a rate read straight off a genetic map would let the path move far more freely than any tuned
 * value has. A map should come with a smoothing factor.
 *
 * @param distances distance covered by the step *into* each position, so element `i` governs the step
 *   from `i - 1` to `i`. Element 0 has no predecessor and is never read by the recursions; it is
 *   filled from element 1 so the array holds no surprises if something else does read it.
 * @param probSwitch probability of at least one switch across [referenceDistance]
 * @param referenceDistance the distance [probSwitch] is quoted over, in the same unit as [distances]
 * @param alternativeCount how many states a switch may go *to*, which is the divisor the recursion
 *   needs: the state count for a haploid recursion, the founder count for the diploid fast path,
 *   where the transition is per haplotype
 * @return the log of staying, and the log of switching to one particular alternative, per position
 */
fun transitionLogsByDistance(
    distances: IntArray,
    probSwitch: Double,
    referenceDistance: Double,
    alternativeCount: Int
): Pair<DoubleArray, DoubleArray> {
    require(distances.isNotEmpty()) { "distances must hold one entry per position" }
    require(probSwitch > 0.0 && probSwitch < 1.0) { "prob-switch must be between 0 and 1 exclusive" }
    require(referenceDistance > 0.0) { "reference distance must be positive" }

    val rate = switchRate(probSwitch, referenceDistance)
    val lnNoSwitch = DoubleArray(distances.size)
    val lnSwitch = DoubleArray(distances.size)

    for (position in distances.indices) {
        // A step of zero -- two records at the same position -- must not be free to switch across,
        // so the distance is floored at one unit.
        val distance = distances[position].coerceAtLeast(1).toDouble()
        val exponent = -rate * distance
        lnNoSwitch[position] = exponent
        lnSwitch[position] = if (alternativeCount > 1) {
            // -expm1(exponent) is 1 - exp(exponent) without the cancellation that subtraction from
            // 1 suffers when the exponent is tiny, which is the normal case here: at 1e-4 per Mb and
            // a 1 kb step the exponent is -1e-7, where the naive form loses most of its precision.
            ln(-expm1(exponent) / (alternativeCount - 1))
        } else {
            // Nowhere to switch to, so switching is impossible rather than free.
            -1.0e6
        }
    }
    if (distances.size > 1) {
        lnNoSwitch[0] = lnNoSwitch[1]
        lnSwitch[0] = lnSwitch[1]
    }
    return Pair(lnNoSwitch, lnSwitch)
}

/**
 * Distance covered by the step into each position, from ascending reference [positions].
 *
 * Element 0 is 0, having no predecessor. A non-ascending input is rejected rather than producing a
 * negative distance, since that would silently become a free switch.
 */
fun stepDistances(positions: IntArray): IntArray {
    val distances = IntArray(positions.size)
    for (index in 1 until positions.size) {
        val step = positions[index] - positions[index - 1]
        require(step >= 0) {
            "positions must ascend, but position $index (${positions[index]}) precedes " +
                    "${positions[index - 1]}"
        }
        distances[index] = step
    }
    return distances
}

/**
 * Expected number of switches across [length] units at [probSwitch] per [referenceDistance] -- the
 * figure that makes a value interpretable, since the parameter itself is quoted over an arbitrary
 * reference distance.
 */
fun expectedSwitches(length: Double, probSwitch: Double, referenceDistance: Double): Double =
    switchRate(probSwitch, referenceDistance) * length

/** Probability of staying on the same founder across [distance]; the complement of a switch. */
fun probabilityOfNoSwitch(distance: Double, probSwitch: Double, referenceDistance: Double): Double =
    exp(-switchRate(probSwitch, referenceDistance) * distance)

/**
 * The switch rate per unit distance, -ln(1 - [probSwitch]) / [referenceDistance]. Computed with ln1p:
 * 1 - p rounds to exactly 1 in double precision once p is below about 1e-16, which would make the rate
 * zero and switching impossible, and that is the range a tuned value can reach.
 */
private fun switchRate(probSwitch: Double, referenceDistance: Double): Double =
    -ln1p(-probSwitch) / referenceDistance
