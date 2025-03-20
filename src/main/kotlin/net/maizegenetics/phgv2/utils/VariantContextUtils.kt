package net.maizegenetics.phgv2.utils

import htsjdk.variant.variantcontext.VariantContext
import org.apache.logging.log4j.LogManager

class VariantContextUtils {
    companion object {
        val myLogger = LogManager.getLogger(VariantContextUtils::class.java)


        /**
         * Function will return -1 if unable to resize the variantcontext due to its type(mostly an INDEL)
         * If the requested position is outside of the current variants coordinates it will return the ASM_Start for + strand and ASM_End for - strand
         */
        fun resizeVariantContext(variant: VariantContext, position: Int, strand : String) : Int {
            //check to see if the variant is either a RefBlock or is a SNP with equal lengths
            return if(isVariantResizable(variant)) {
                when {
                    position < variant.start -> variant.getAttributeAsInt("ASM_Start",variant.start)
                    position > variant.end -> variant.getAttributeAsInt("ASM_End",variant.end)
                    strand == "+" -> {
                        val offset = position - variant.start
                        variant.getAttributeAsInt("ASM_Start",variant.start) + offset
                    }
                    strand == "-" -> {
                        val offset = position - variant.start
                        variant.getAttributeAsInt("ASM_Start",variant.end) - offset
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
        fun isVariantResizable(variant: VariantContext) : Boolean {
            return when {
                variant.getReference().baseString.length == 1 && variant.end - variant.start > 0 && variant.type == VariantContext.Type.SYMBOLIC -> true //refBlock
                variant.reference.baseString.length == variant.getAlternateAllele(0).baseString.length -> true //This covers both SNPs and multiallelic polymorphisms
                else -> false
            }
        }

        /**
         * Function to see if the BED region is fully contained within a VariantContext
         * Indels are left-justified
         * Bed:       |---|
         * Var: |--------------|
         */
        fun bedRegionContainedInVariant(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
            val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
            return variant.contig == region.first.contig && variant.start <= region.first.position && end >= region.second.position
        }

        /**
         * Function to see if the VariantContext is fully contained within a BED region
         * Indels are left-justified
         * Bed: |--------------|
         * Var:       |---|
         */
        fun variantFullyContained(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
            val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
            return variant.contig == region.first.contig && variant.start >= region.first.position && end <= region.second.position
        }

        /**
         * Function to see if the start of the variant is partially contained in the BED region
         * Indels are left-justified
         * Bed: |--------------|
         * Var:              |---|
         */
        fun variantPartiallyContainedStart(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
            val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
            return variant.contig == region.first.contig &&
                    variant.start >= region.first.position &&
                    variant.start <= region.second.position &&
                    end > region.second.position
        }

        /**
         * Function to see if the end of the variant is partially contained in the BED region
         * Indels are left-justified
         * Bed:      |--------------|
         * Var:    |---|
         */
        fun variantPartiallyContainedEnd(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
            val end = if (variant.type == VariantContext.Type.SYMBOLIC) variant.end else variant.start
            return variant.contig == region.first.contig &&
                    end <= region.second.position &&
                    end >= region.first.position &&
                    variant.start < region.first.position
        }

        /**
         * Function to see if the variant is after the bed region
         * Bed: |--------------|
         * Var:                   |---|
         */
        fun variantAfterRegion(region: Pair<Position,Position>, variant: VariantContext) : Boolean {
            return variant.contig == region.first.contig && variant.start > region.second.position
        }


    }
}