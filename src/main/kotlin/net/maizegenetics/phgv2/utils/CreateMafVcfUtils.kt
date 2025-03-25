package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeaderLine
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.cli.HVCFRecordMetadata
import org.apache.logging.log4j.LogManager


class CreateMafVcfUtils {
    companion object {
        val myLogger = LogManager.getLogger(CreateMafVcfUtils::class.java)


        /**
         * Function to resize the Position range based on the new position.  If isFirst is true then it will resize the start position, otherwise it will resize the end position
         */
        fun resizePositionRange(positionRange: Pair<Position,Position>, newPosition : Int, isFirst: Boolean) : Pair<Position,Position> {
            return if(isFirst) {
                //Slide the start position to the new position
                Pair(Position(positionRange.first.contig,newPosition),positionRange.second)
            } else {
                //Slide the end position to the new position
                Pair(positionRange.first,Position(positionRange.second.contig,newPosition))
            }
        }


        /**
         * This function will bulk load sequences in from the AGC record and then will associate the returned sequences
         * with the metadata record which contains the coordiantes for the query and will add in the asmSeq.
         */
        fun addSequencesToMetaData(dbPath: String, metadata: List<HVCFRecordMetadata>, condaEnvPrefix:String="") : List<HVCFRecordMetadata> {
            //get out the assembly coordinates and build them into the regions
            val metaDataToRangeLookup = metadata.map {
                val queries = mutableListOf<String>()
                val displayNames = mutableListOf<String>()

                for(range in it.asmRegions) {
                    if(range.first.position-1 > range.second.position-1) {
                        queries.add("${range.first.contig}@${it.sampleName}:${range.second.position-1}-${range.first.position-1}")
                        displayNames.add("${range.first.contig}:${range.second.position-1}-${range.first.position-1}")
                    }
                    else {
                        queries.add("${range.first.contig}@${it.sampleName}:${range.first.position - 1}-${range.second.position - 1}")
                        displayNames.add("${range.first.contig}:${range.first.position-1}-${range.second.position-1}")
                    }
                }

                Triple(it,queries, displayNames)
            }

            val ranges = metaDataToRangeLookup.flatMap { it.second }

            val seqs = retrieveAgcContigs(dbPath,ranges,condaEnvPrefix)

            return metaDataToRangeLookup.map { it.first.copy(asmSeq = buildSeq(seqs,it.third,it.first)) } //This is a useful way to keep things immutable
        }


        /**
         * Function to build the haplotype sequence based on the list of display regions and the given haplotype sequence object.
         * The sequence is already extracted out of AGC and stored in the seqs map.
         * The seqs map is keyed by a Pair of (sampleName, displayRegion) and the value is the NucSeq object.
         * The Pair may look something like ("B97", "1:1-1000") or ("B97", "chr1")
         */
        fun buildSeq(seqs: Map<Pair<String,String>, NucSeq>, displayRegions : List<String>, hvcfRecordMetadata: HVCFRecordMetadata) : String {
            val hapSeqRegions = hvcfRecordMetadata.asmRegions

            return displayRegions.mapIndexed{ idx, currentDisplayRegion ->
                val currentHapSeqRegion = hapSeqRegions[idx]

                val seq = seqs[Pair(hvcfRecordMetadata.sampleName,currentDisplayRegion)]!!

                //Means it is the first region
                if(currentHapSeqRegion.first.position > currentHapSeqRegion.second.position) {
                    //Means it needs to be reverse complemented
                    seq.reverse_complement().seq()
                }
                else {
                    seq.seq()
                }
            }.joinToString()
        }

        /**
         * Simple function to convert the all the HVCFRecordMetadata records into VariantContext records
         */
        fun convertMetaDataToHVCFContexts(metaData: List<HVCFRecordMetadata>, asmHeaders: MutableMap<String, VCFHeaderLine>, dbPath:String): List<VariantContext> {
            return metaData.map { convertMetaDataRecordToHVCF(it, asmHeaders, dbPath) }
        }

        /**
         * Simple function to convert a single metadata record into a VariantContext.
         * This will also create the ALT tag and add it to the asmHeaders object for use during export.
         */
        private fun convertMetaDataRecordToHVCF(metaDataRecord: HVCFRecordMetadata, asmHeaders: MutableMap<String, VCFHeaderLine>, dbPath: String): VariantContext {
            val assemblyHaplotypeSeq:String = metaDataRecord.asmSeq
            //md5 hash the assembly sequence
            val assemblyHaplotypeHash = getChecksumForString(assemblyHaplotypeSeq)
            check(metaDataRecord.refSeq.isNotEmpty()) { "Reference sequence is empty" }
            //md5 has the refSequence
            val refSeqHash = getChecksumForString(metaDataRecord.refSeq)

            //create the asmHeader lines
            if(!asmHeaders.containsKey(assemblyHaplotypeHash)) {
                asmHeaders[assemblyHaplotypeHash] =
                    VCFAltHeaderLine(
                        "<ID=${assemblyHaplotypeHash}, Description=\"haplotype data for line: ${metaDataRecord.sampleName}\">," +
                                "Source=\"${dbPath}/assemblies.agc\",SampleName=\"${metaDataRecord.sampleName}\"," +
                                "Regions=\"${metaDataRecord.asmRegions.map { "${it.first.contig}:${it.first.position}-${it.second.position}" }.joinToString(",")}\"," +
                                "Checksum=\"${assemblyHaplotypeHash}\",RefChecksum=\"${refSeqHash}\",RefRange=\"${metaDataRecord.refContig}:${metaDataRecord.refStart}-${metaDataRecord.refEnd}\">",
                        VCFHeaderVersion.VCF4_2
                    )
            } else {
                myLogger.info("convertMetaDataRecordToHVCF: asmHeaders already contains key ${assemblyHaplotypeHash}")
            }


            //build a variant context of the HVCF with the hashes
            return createHVCFRecord(metaDataRecord.sampleName, Position(metaDataRecord.refContig,metaDataRecord.refStart),
                Position(metaDataRecord.refContig, metaDataRecord.refEnd ),
                Pair(metaDataRecord.refSeq[0].toString(), assemblyHaplotypeHash))
        }

        /**
         * Strand aware function to merge together consecutive assembly regions.  This is done to reduce the number of entries in the hvcf alt header.
         */
        fun mergeConsecutiveRegions(variants: List<Pair<Position,Position>>) : List<Pair<Position,Position>> {
            val mergedConsecutiveVariants = mutableListOf<Pair<Position,Position>>()
            var currentStart = variants.first().first
            var currentEnd = variants.first().second
            for(i in 1 until variants.size) {
                val nextStart = variants[i].first
                val nextEnd = variants[i].second

                //Check to see if the next region is only 1 bp long.  If so we need to check both normal and inverted boundaries and extend if it makes sense, if not reset the currentStart and currentEnd
                if(nextStart.position == nextEnd.position) {
                    //Check to see if the next region is on the + strand
                    //Using nextStart here as it equals nextEnd
                    if(nextStart.position == currentEnd.position + 1 || nextStart.position == currentEnd.position -1) {
                        currentEnd = nextStart
                    }
                    else {
                        mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
                        currentStart = nextStart
                        currentEnd = nextEnd
                    }
                }
                else if(currentEnd < currentStart && nextEnd < nextStart) {
                    //This is the case where we have a variant that is on the - strand
                    //We need to check if the next variant is consecutive
                    if(nextStart.position == (currentEnd.position - 1)) {
                        currentEnd = nextEnd
                    }
                    else {
                        mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
                        currentStart = nextStart
                        currentEnd = nextEnd
                    }
                }
                else if(nextStart.position == (currentEnd.position + 1)) {
                    currentEnd = nextEnd
                }
                else {
                    mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
                    currentStart = nextStart
                    currentEnd = nextEnd
                }
            }
            mergedConsecutiveVariants.add(Pair(currentStart,currentEnd))
            return mergedConsecutiveVariants
        }
    }
}