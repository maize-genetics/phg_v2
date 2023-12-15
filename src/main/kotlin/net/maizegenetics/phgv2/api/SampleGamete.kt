package net.maizegenetics.phgv2.api

/**
 * A SampleGamete is a combination of a sample name and a gamete id.
 */
data class SampleGamete(val name: String, val gameteId: Int = 0) {
}