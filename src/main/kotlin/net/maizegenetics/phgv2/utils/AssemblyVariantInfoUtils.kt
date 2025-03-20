package net.maizegenetics.phgv2.utils
import biokotlin.genome.AssemblyVariantInfo
import org.apache.logging.log4j.LogManager

class AssemblyVariantInfoUtils {
    companion object {
        val myLogger = LogManager.getLogger(AssemblyVariantInfoUtils::class.java)

        /**
         * Function to see if the BED region is fully contained within a VariantInfo
         * Indels are left-justified
         * Bed:       |---|
         * Var: |--------------|
         */
        fun bedRegionContainedInVariantInfo(region: Pair<Position,Position>, variant: AssemblyVariantInfo) : Boolean {
            return variant.chr == region.first.contig && variant.startPos <= region.first.position && variant.endPos >= region.second.position
        }



        /**
         * Function to see if the VariantContext is fully contained within a BED region
         * Indels are left-justified
         * Bed: |--------------|
         * Var:       |---|
         */
        fun variantInfoFullyContained(region: Pair<Position,Position>, variant: AssemblyVariantInfo) : Boolean {
            return variant.chr == region.first.contig && variant.startPos >= region.first.position && variant.endPos <= region.second.position
        }


        /**
         * Function to see if the start of the variant is partially contained in the BED region
         * Indels are left-justified
         * Bed: |--------------|
         * Var:              |---|
         */
        fun variantInfoPartiallyContainedStart(region: Pair<Position,Position>, variant: AssemblyVariantInfo) : Boolean {
            return variant.chr == region.first.contig &&
                    variant.startPos >= region.first.position &&
                    variant.startPos <= region.second.position &&
                    variant.endPos > region.second.position
        }



        /**
         * Function to see if the end of the variant is partially contained in the BED region
         * Indels are left-justified
         * Bed:      |--------------|
         * Var:    |---|
         */
        fun variantInfoPartiallyContainedEnd(region: Pair<Position,Position>, variant: AssemblyVariantInfo) : Boolean {
            return variant.chr == region.first.contig &&
                    variant.endPos <= region.second.position &&
                    variant.endPos >= region.first.position &&
                    variant.startPos < region.first.position
        }



        /**
         * Function to see if the variant is after the bed region
         * Bed: |--------------|
         * Var:                   |---|
         */
        fun variantInfoAfterRegion(region: Pair<Position,Position>, variant: AssemblyVariantInfo) : Boolean {
            return variant.chr == region.first.contig && variant.startPos > region.second.position
        }
    }
}