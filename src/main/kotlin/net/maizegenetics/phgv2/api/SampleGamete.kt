package net.maizegenetics.phgv2.api

/**
 * A SampleGamete is a combination of a sample name and a gamete id.
 */
data class SampleGamete(val name: String, val gameteId: Int = 0) : Comparable<SampleGamete> {
    override fun compareTo(other: SampleGamete): Int {
        val namesEqual = name.compareTo(other.name)
        return if (namesEqual == 0) gameteId.compareTo(other.gameteId) else namesEqual
    }

    override fun equals(other: Any?): Boolean {
        return if (other is SampleGamete) {
            name == other.name && gameteId == other.gameteId
        } else false
    }

    override fun toString(): String {
        return "$name:$gameteId"
    }

    override fun hashCode(): Int {
        return name.hashCode() + gameteId.hashCode()
    }
}