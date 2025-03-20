package net.maizegenetics.phgv2.utils

import htsjdk.variant.variantcontext.VariantContext
import org.apache.logging.log4j.LogManager

class VariantContextUtils {
    companion object {
        val myLogger = LogManager.getLogger(VariantContextUtils::class.java)

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