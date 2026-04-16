package net.maizegenetics.phgv2.utils


import biokotlin.genome.SeqRangeSort
import biokotlin.seq.NucSeq
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.HVCFRecordMetadata
import org.apache.logging.log4j.LogManager
import java.io.File


data class HvcfVariant(val refRange: ReferenceRange,val sampleName: String, val hapId: String)


class VCFConversionUtils {
    companion object {
        val myLogger = LogManager.getLogger(VCFConversionUtils::class.java)

        /**
         * Function to convert a GVCF file into an HCVF file
         */
        fun convertGVCFToHVCF(dbPath: String, sampleName: String, bedRanges : List<Pair<Position,Position>>, gvcfVariants: List<VariantContext>,
                              refGenomeSequence : Map<String, NucSeq>, agcArchiveName: String, asmHeaders: MutableMap<String, VCFHeaderLine>, condaEnvPrefix:String = "", numRangesPerAgcQuery: Int = -1) : List<VariantContext> {
            // group the gvcfVariants by contig
            val gvcfVariantsByContig = gvcfVariants.groupBy { it.contig }

            val bedRegionsByContig = bedRanges.groupBy { it.first.contig }


            myLogger.info("in convertGVCFToHVCF: sort and call converGVCFToHVCFForChrom")
            return gvcfVariantsByContig.keys
                .sortedWith(compareBy(SeqRangeSort.alphaThenNumberSort){ name:String -> name}) //Need to do a sort here as we need to make sure we process the chromosomes in
                .filter { bedRegionsByContig.containsKey(it) }
                .flatMap { convertGVCFToHVCFForChrom(dbPath, sampleName, bedRegionsByContig[it]!!, refGenomeSequence, agcArchiveName, gvcfVariantsByContig[it]!!, asmHeaders,condaEnvPrefix, numRangesPerAgcQuery) }
        }

        private fun convertGVCFToHVCFForChrom(dbPath: String, sampleName: String, bedRanges: List<Pair<Position,Position>>, refGenomeSequence: Map<String, NucSeq>, agcArchiveName: String, variantContexts: List<VariantContext>, asmHeaders: MutableMap<String, VCFHeaderLine>, condaEnvPrefix:String = "", numRangesPerAgcQuery: Int = -1 ) : List<VariantContext> {

            /**
             * Loop through the bed file
             * Loop through the gvcf records as well
             *
             * We need to determine if our BED region overlaps with the gvcf record
             * To do this we need to collect the gvcf records into their corresponding bed Regions
             * Then from those collected regions, we take the first and last ones and resize based on the bed regions to get the asm_Start and asm_End
             * Then extract the sequence out of the AGC archive and md5 hash it
             * Then call the createHVCFRecord with this information
             */
            myLogger.info("in convertGVCFToHVCFForChrom: bedRanges.size = ${bedRanges.size}")
            val outputVariantMetadata = mutableListOf<HVCFRecordMetadata>()
            var currentVariantIdx = 0
            for(region in bedRanges) {
                val regionStart = region.first.position
                val regionEnd = region.second.position
                val regionChrom = region.first.contig
                val tempVariants = mutableListOf<VariantContext>()

                check(regionChrom in refGenomeSequence.keys) { "Chromosome $regionChrom not found in reference" }

                //Need to subtract here as the Biokotlin NucSeq is 0 based
                val refRangeSeq = refGenomeSequence[regionChrom]!![regionStart-1..regionEnd-1]

                while (currentVariantIdx < variantContexts.size) {
                    val currentVariant = variantContexts[currentVariantIdx]

                    //check different cases for the variant
                    //If variant is fully contained in Bed region add to temp list and increment currentVariantIdx
                    //If variant is partially contained in Bed region add to temp list do not increment as we need to see if the next bed also overlaps
                    //If variant is not contained in Bed region, skip and do not increment as we need to see if the next bed overlaps
                    if(VariantContextUtils.bedRegionContainedInVariant(region, currentVariant)) {
                        outputVariantMetadata.add(
                            convertGVCFRecordsToHVCFMetaData(
                                sampleName,
                                region,
                                refRangeSeq,
                                listOf(currentVariant)
                            )
                        )
                        tempVariants.clear()
                        break
                    }
                    if(VariantContextUtils.variantFullyContained(region, currentVariant)) {
                        //This is the case where the variant is completely contained within the region
                        tempVariants.add(currentVariant)
                        currentVariantIdx++
                    }
                    else if(VariantContextUtils.variantPartiallyContainedStart(region,currentVariant)) {
                        tempVariants.add(currentVariant)
                        break
                    }
                    else if(VariantContextUtils.variantPartiallyContainedEnd(region, currentVariant)) {
                        tempVariants.add(currentVariant)
                        currentVariantIdx++
                    }
                    else if(VariantContextUtils.variantAfterRegion(region, currentVariant)) {
                        //write out what is in tempVariants
                        if(tempVariants.isNotEmpty()) {
                            outputVariantMetadata.add(
                                convertGVCFRecordsToHVCFMetaData(sampleName,
                                    region,
                                    refRangeSeq,
                                    tempVariants
                                )
                            )

                            tempVariants.clear()
                        }
                        //move up Bed region
                        break
                    }
                    else { //this is the case if the Variant is behind the BED region
                        //move up Variant
                        currentVariantIdx++
                    }
                }

                if(tempVariants.isNotEmpty()) {
                    outputVariantMetadata.add(convertGVCFRecordsToHVCFMetaData(
                        sampleName,
                        region,
                        refRangeSeq,
                        tempVariants
                    ))
                    tempVariants.clear()
                }
            }

            val metaDataWithSequence = CreateMafVcfUtils.addSequencesToMetaData(dbPath, outputVariantMetadata, condaEnvPrefix, numRangesPerAgcQuery)
            val outputVariants = CreateMafVcfUtils.convertMetaDataToHVCFContexts(metaDataWithSequence, asmHeaders, dbPath)

            return outputVariants
        }

        /**
         * Function to extract all the needed information out of the ASM gVCF record and put them into HVCFRecordMetadata objects
         * This will first try to resize the positions based on the ref start position and then will extract out all the other information.
         */
        private fun convertGVCFRecordsToHVCFMetaData(sampleName: String, region: Pair<Position,Position>, refRangeSeq: NucSeq, variants: List<VariantContext> ) : HVCFRecordMetadata {
            //Take the first and the last variantContext
            val firstVariant = variants.first()
            val lastVariant = variants.last()

            //val check strandedness of the variants
            val firstStrand = firstVariant.getAttributeAsString("ASM_Strand","+")

            val lastStrand = lastVariant.getAttributeAsString("ASM_Strand","+")
            //Resize the first and last variantContext ASM start and end based on the regions
            var newASMStart = VariantContextUtils.resizeVariantContext(firstVariant, region.first.position, firstStrand)
            if(newASMStart == -1) {
                newASMStart = if(firstStrand == "+") firstVariant.getAttributeAsInt("ASM_Start",region.first.position)
                else firstVariant.getAttributeAsInt("ASM_End",region.first.position)
            }

            var newASMEnd = VariantContextUtils.resizeVariantContext(lastVariant, region.second.position, lastStrand)
            if(newASMEnd == -1) {
                newASMEnd = if(lastStrand == "+") lastVariant.getAttributeAsInt("ASM_End",region.second.position)
                else lastVariant.getAttributeAsInt("ASM_Start",region.second.position)
            }

            val regions = buildNewAssemblyRegions(newASMStart,newASMEnd,variants)


            return HVCFRecordMetadata(sampleName=sampleName, refSeq = refRangeSeq.toString(), asmSeq = "",
                refContig = region.first.contig, refStart = region.first.position, refEnd = region.second.position,
                regions)

        }





        /**
         * Function to build the new assembly region coordinates based on the new Start and end and the list of VariantContexts
         * Any consecutive regions should be merged together so we do not make the eventual string too long
         * The output will be a List<Pair<Position,Position>> which will be the new coordinates for the all the  assembly regions
         */
        fun buildNewAssemblyRegions(newStart: Int, newEnd: Int, variants: List<VariantContext>) : List<Pair<Position,Position>> {
            val variantsConverted = variants.map { convertVariantContextToPositionRange(it) }

            //resize the first and last position based on the strand
            val resizedFirst = CreateMafVcfUtils.resizePositionRange(variantsConverted.first(),newStart,true)
            val resizedLast = CreateMafVcfUtils.resizePositionRange(variantsConverted.last(),newEnd,false)

            //merge the first and last with the rest of the variants
            val mergedVariants = mutableListOf<Pair<Position,Position>>()
            if(variantsConverted.size == 1) {
                val resizedFirstAndLast = CreateMafVcfUtils.resizePositionRange(CreateMafVcfUtils.resizePositionRange(variantsConverted.first(), newStart, true), newEnd, false)
                mergedVariants.add(resizedFirstAndLast)
            } else {
                mergedVariants.add(resizedFirst) // add the first variant
                if (variantsConverted.size > 2) { // add any variants in the middle
                    mergedVariants.addAll(variantsConverted.subList(1, variantsConverted.size - 1))
                }
                mergedVariants.add(resizedLast) // add the last variant
            }

            //merge the consecutive regions
            val mergedConsecutiveVariants =CreateMafVcfUtils.mergeConsecutiveRegions(mergedVariants)

            return mergedConsecutiveVariants
        }

        /**
         * Function to convert a VariantContext into a Pair<Position,Position> which will be the assembly starts and ends of the variantContext
         */
        fun convertVariantContextToPositionRange(variant: VariantContext) : Pair<Position,Position> {
            //get out the assembly coords
            val contig = variant.getAttributeAsString("ASM_Chr","")
            val start = variant.getAttributeAsInt("ASM_Start",variant.start)
            val end = variant.getAttributeAsInt("ASM_End",variant.end)
            return Pair(Position(contig, start), Position(contig, end))
        }


        //The Map<Pair<RefRange,HapId>,List<HvcfVariant>>
        // the list does not matter as the sequence is the same but we should collect it in case
        //We retain all of them as we will need those boundaries later and it is easier to just keep them together
        fun createASMHapIdMap(dbPath: String): Map<Pair<ReferenceRange,String>, List<HvcfVariant>> {
            //check to see if there are hvcfFiles in the dbPath

            //If so we need to loop through them and extract out the HapIds and pair them with the refRanges.
            //TODO Eventually we can replace this with direct tileDB calls
            return File("${dbPath}/hvcf_files").walk().filter { !it.isHidden && !it.isDirectory }
                .filter {
                    it.name.endsWith(".h.vcf.gz") || it.name.endsWith(".h.vcf") ||
                            it.name.endsWith(".hvcf.gz") || it.name.endsWith(".hvcf")
                }
                .flatMap { hvcfFile ->
                    extractASMHapIds(hvcfFile)
                }
                .groupBy { Pair(it.refRange, it.hapId) } // We then can collect across all samples
        }

        fun createASMNameAndRefRangeMap(dbPath: String): Map<Pair<ReferenceRange,String>, List<HvcfVariant>> {
            //TODO Eventually we can replace this with direct tileDB calls
            return File("${dbPath}/hvcf_files").walk().filter { !it.isHidden && !it.isDirectory }
                .filter {
                    it.name.endsWith(".h.vcf.gz") || it.name.endsWith(".h.vcf") ||
                            it.name.endsWith(".hvcf.gz") || it.name.endsWith(".hvcf")
                }
                .flatMap { hvcfFile ->
                    extractASMHapIds(hvcfFile)
                }
                .groupBy { Pair(it.refRange, it.sampleName) } // We then can collect across all samples
        }

        /**
         * Function to pull out the assembly hvcf information so we can look it up later.
         * Here we can use the first alt allele as it is guaranteed to be haploid because it is coming from tiledb.
         */
        fun extractASMHapIds(file: File): List<HvcfVariant> {
            val hvcfReader = VCFFileReader(file, false)
            val sampleName = hvcfReader.header.sampleNamesInOrder.first()
            val hvcfVariants = hvcfReader.mapNotNull { context ->
                if (context.alternateAlleles.isEmpty() || context.alternateAlleles.any { it.isNoCall }) {
                    null
                } else {
                    //Here we can get the first alt allele as the ASM hvcfs are guaranteed to be single sample haploid
                    val hapId = context.getAlternateAllele(0).displayString.removeSurrounding("<", ">")
                    val refRange = ReferenceRange(context.contig, context.start, context.end)
                    HvcfVariant(refRange, sampleName,hapId)
                }
            }

            return hvcfVariants
        }
    }
}