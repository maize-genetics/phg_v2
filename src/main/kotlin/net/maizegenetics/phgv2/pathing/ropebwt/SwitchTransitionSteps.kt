package net.maizegenetics.phgv2.pathing.ropebwt

import kotlin.math.exp

/**
 * Forward-backward transition steps for the switch model [ImputePathFromPs4g] uses: each haplotype
 * stays on its parent or switches to any one of the other parents with equal probability, at
 * probabilities that vary by position (see [transitionLogsByDistance]).
 *
 * The step is applied without materialising a matrix. Stay is `s` and switching to one particular
 * other parent is `w`, so for a vector `v` over the parents with total `S`,
 *
 *     out[j] = s v[j] + w (S - v[j])
 *
 * which is O(nParents) per step. The matrix is symmetric, so the forward and backward steps are the
 * same.
 *
 * For a [diploid] path the two haplotypes switch independently -- an inbreeding coefficient of 0, as
 * in [ViterbiHMM] -- so the transition over ordered pairs is the Kronecker product of two haploid
 * transitions, `T((a, b) -> (c, d)) = t(a, c) t(b, d)`. It is applied as the haploid step along the
 * second founder for each first founder, then along the first founder for each second, which is
 * O(nParents^2) per step rather than O(nParents^4) for the full matrix.
 *
 * @param lnNoSwitch per-position log probability that one haplotype stays on its parent.
 * @param lnSwitch per-position log probability that one haplotype switches to one particular other parent.
 */
class SwitchTransitionStep(
    private val nParents: Int,
    private val diploid: Boolean,
    lnNoSwitch: DoubleArray,
    lnSwitch: DoubleArray
) : TransitionStep {

    private val stay = DoubleArray(lnNoSwitch.size) { exp(lnNoSwitch[it]) }
    private val switch = DoubleArray(lnSwitch.size) { exp(lnSwitch[it]) }
    private val buffer = DoubleArray(if (diploid) nParents * nParents else 0)

    init {
        require(lnNoSwitch.size == lnSwitch.size) { "lnNoSwitch and lnSwitch must have the same length" }
    }

    override fun apply(position: Int, vector: DoubleArray, out: DoubleArray, forward: Boolean) {
        val s = stay[position]
        val w = switch[position]
        if (!diploid) {
            haploidStep(s, w, vector, 0, 1, out)
            return
        }
        // along the second founder: buffer[a, d] = sum over b of vector[a, b] t(b, d)
        for (first in 0 until nParents) haploidStep(s, w, vector, first * nParents, 1, buffer)
        // along the first founder: out[c, d] = sum over a of buffer[a, d] t(a, c)
        for (second in 0 until nParents) haploidStep(s, w, buffer, second, nParents, out)
    }

    /**
     * The haploid step over the nParents entries of [vector] starting at [offset] and spaced
     * [stride] apart, written to the same entries of [out].
     */
    private fun haploidStep(s: Double, w: Double, vector: DoubleArray, offset: Int, stride: Int, out: DoubleArray) {
        var total = 0.0
        for (i in 0 until nParents) total += vector[offset + i * stride]
        for (i in 0 until nParents) {
            val index = offset + i * stride
            out[index] = s * vector[index] + w * (total - vector[index])
        }
    }
}
