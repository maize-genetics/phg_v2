package net.maizegenetics.phgv2.utils
import biokotlin.genome.AssemblyVariantInfo
import org.apache.logging.log4j.LogManager

class AssemblyVariantInfoUtils {
    companion object {
        val myLogger = LogManager.getLogger(AssemblyVariantInfoUtils::class.java)


        /**
         * Function will return -1 if unable to resize the variantcontext due to its type(mostly an INDEL)
         * If the requested position is outside of the current variants coordinates it will return the ASM_Start for + strand and ASM_End for - strand
         */
        fun resizeVariantInfo(variant: AssemblyVariantInfo, position: Int, strand : String) : Int {
            //check to see if the variant is either a RefBlock or is a SNP with equal lengths
            return if(isVariantInfoResizable(variant)) {
                when {
                    position < variant.startPos -> variant.asmStart
                    position > variant.endPos -> variant.asmEnd
                    strand == "+" -> {
                        val offset = position - variant.startPos
                        variant.asmStart + offset
                    }
                    strand == "-" -> {
                        val offset = position - variant.startPos
                        variant.asmStart - offset
                    }
                    else -> -1
                }
            } else {
                -1
            }
        }



        /**
         * Function to check if a variant is resizable.  Only RefBlocks and SNPs are resizable
         */
        fun isVariantInfoResizable(variant: AssemblyVariantInfo) : Boolean {
            return when {
                variant.refAllele.length == 1 && variant.endPos - variant.startPos > 0 -> true //refBlock
                variant.refAllele.length == variant.altAllele.length -> true //This covers both SNPs and multiallelic polymorphisms
                else -> false
            }
        }

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