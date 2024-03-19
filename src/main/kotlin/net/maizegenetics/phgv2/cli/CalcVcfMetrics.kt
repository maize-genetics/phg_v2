package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.google.common.collect.DiscreteDomain
import com.google.common.collect.Range
import com.google.common.collect.TreeRangeSet
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.Exception
import kotlin.math.abs

// values calculated for GVCF files
// single chrom, or "ALL" for all contigs
data class GVCFStats(
    val sample: String,
    val chrom: String,
    val refLength: Int,
    val numSNPs: Int,
    val numIns: Int,
    val numDel: Int,
    val numNs: Int,
    val numBasesInserted: Int,
    val numBasesDeleted: Int,
    val percentIdentityWithRef: Double,
    val percentMappedToRef: Double,
    val meanInsertionSize: Double,
    val medianInsertionSize: Double?,
    val largestInsertion: Int?,
    val meanDeletionSize: Double,
    val medianDeletionSize: Double?,
    val largestDeletion: Int?
) {
    override fun toString(): String {
        return "$sample\t$chrom\t$refLength\t$numSNPs\t$numIns\t$numDel\t$numNs\t$numBasesInserted\t$numBasesDeleted\t$percentIdentityWithRef\t" +
                "$percentMappedToRef\t$meanInsertionSize\t$medianInsertionSize\t$largestInsertion\t$meanDeletionSize\t$medianDeletionSize\t" +
                "$largestDeletion"
    }

}

data class HVCFStatsByRange(
    val sampleID: String,
    val sampleLength: Int,
    val identicalToRef: Boolean
)
data class HVCFStats(
    val sample: String,
    val chrom: String,
    val refRangesWithHap: Int,
    val hapsIdenticalToRef: Int
)

data class RefRangeInfo (
    val id: String,
    val chrom: String,
    val start: Int,
    val end: Int
        )

data class BaseCounter(
    val rangesMapped: TreeRangeSet<Int>, // ranges covered by some gvcf record, including missing
    val ins: MutableList<Int>, // sizes of insertions
    val del: MutableList<Int>, // sizes of deletions
    var snps: Int, // number of SNPs
    var ns: Int, // number of Ns - includes indels
    var ref: Int // number of reference bases
)



/**
 * A [CliktCommand] class for calculating metrics on g.vcf files to be used in quality control
 */
class CalcVcfMetrics: CliktCommand(help="Calculate quality control metrics on g.vcf files") {

    private val myLogger = LogManager.getLogger(LoadVcf::class.java)

    val vcfDir by option("--vcf-dir", help = "Full path to VCF file directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--vcf-dir must not be blank"
            }
        }

    val gvcfOutFile by option( "--gvcf-output", help = "Name for the GVCF output .tsv file. If blank, stats will not be calculated for g.vcf files.")
        .default("")

    override fun run() {
        calculateVcfMetrics(vcfDir, gvcfOutFile)
    }


    fun calculateVcfMetrics(vcfDir: String, gvcfOutFile: String = "") {
        // get vcf files from directory
        val fileList = File(vcfDir).listFiles()
        val gvcfFileList = fileList.filter{it.isFile && (it.name.endsWith("g.vcf.gz") || it.name.endsWith(".gvcf.gz"))}
        val hvcfFileList = fileList.filter{it.isFile && (it.name.endsWith("h.vcf.gz") || it.name.endsWith(".hvcf.gz"))}

        // map file name to gvcf stats
        val gvcfStats = if(gvcfOutFile != "") {
            if(gvcfFileList.isEmpty()) {
                myLogger.warn("No files ending in g.vcf.gz found in $vcfDir.  " +
                        "Note that both the bgzipped and indexed files must exist in the specified folder. \n" +
                        "Please check the folder and try again.")
                throw IllegalArgumentException("CalcVcfMetrics: No files ending in g.vcf.gz found in $vcfDir.  " +
                        "Note that both the bgzipped and indexed files must exist in the specified folder. \n" +
                        "Please check the folder and try again.")
            } else {
                gvcfFileList.associate{Pair(it.name.removeSuffix(".gvcf.gz").removeSuffix(".g.vcf.gz"), getGVCFStats(it))}
            }
        } else { null }

        // map file name to hvcf stats, both types
        val hvcfStatsByChrom =  if(hvcfFileList.isEmpty()) {
            myLogger.warn("No files ending in h.vcf.gz found in $vcfDir.  " +
                        "Note that both the bgzipped and indexed files must exist in the specified folder. \n" +
                        "Please check the folder and try again.")
            mutableMapOf()
        } else {

            hvcfFileList.associate{
                val fileName = it.name.removeSuffix(".hvcf.gz").removeSuffix(".h.vcf.gz")

                val statsByRange = getHVCFStatsByRange(it)
                Pair(fileName, getHVCFStats(statsByRange.first, statsByRange.second))

            }
        }


        if(gvcfStats != null) {
            writeGVCFFile(gvcfOutFile, gvcfStats, hvcfStatsByChrom)
        }


    }

    /**
     * Takes a list of
     */
    fun writeGVCFFile(outFileName: String, gvcfStatsMap: Map<String, List<GVCFStats>?>, hvcfStatsMap: Map<String, List<HVCFStats>>) {

        val includeHVCFStats = hvcfStatsMap.size > 0

        File(outFileName).bufferedWriter().use {writer ->
            val samples = gvcfStatsMap.keys

            // write header
            writer.write("taxa\tchrom\trefLength\tnumSNPs\tnumIns\tnumDel\tnumNs\tnumBasesInserted\tnumBasesDeleted\t" +
                    "percentIdentityWithRef\tpercentMappedToRef\tmeanInsertionSize\tmedianInsertionSize\tlargestInsertion\t" +
                    "meanDeletionSize\tmedianDeletionSize\tlargestDeletion")

            if(includeHVCFStats) {
                writer.write("\trefRangesWithHaplotype\thaplotypesIdenticalToRef\n")
            } else { writer.write("\n") }

            samples.forEach {sample ->

                val gvcfStats = gvcfStatsMap[sample]!!
                val hvcfStats = hvcfStatsMap[sample]

                gvcfStats.forEach{chrom ->
                    writer.write(chrom.toString())

                    if(includeHVCFStats){
                        val hvcfChrom = hvcfStats?.filter{it.chrom == chrom.chrom}?.get(0)
                        writer.write("\t${hvcfChrom?.refRangesWithHap}\t${hvcfChrom?.hapsIdenticalToRef}\n")
                    } else {
                        writer.write("\n")
                    }
                }
            }
        }
    }


    fun getGVCFStats(vcfFile: File): List<GVCFStats>? {

        VCFFileReader(vcfFile, false).use { reader ->
                val header = reader.fileHeader
                val records = reader.toList()

            return getGVCFStats(header, records)
        }
    }

    fun getGVCFStats(header: VCFHeader, records: List<VariantContext>): List<GVCFStats>? {

            // gvcf files should be single-sample
            if(header.genotypeSamples.size != 1) {
                myLogger.warn("${header.genotypeSamples.size} samples detected. GVCF must be single sample.")
                return null
            }

            val sampleName = header.genotypeSamples[0]
            myLogger.info("Processing sample $sampleName")

            // get all reference chromosomes and their respective lengths
            // includes scaffolds and contigs, where present
            val refChrLengths = header.contigLines
                .filter{it.genericFields["ID"] != null && it.genericFields["length"] != null}
                .associate { Pair(it.genericFields["ID"], it.genericFields["length"]!!.toInt())}


            // keep track of indels, SNPs, Ns, by contig name
            val baseCounterMap = mutableMapOf<String, BaseCounter>()

            // iterate through gcvf records
            val iterator = records.iterator()

            while(iterator.hasNext()) {
                val record = iterator.next()

                // chromosome
                val refChrom = record.contig

                // asm start and end positions, mostly so that we don't have to type the long call multiple times
                val asmEnd = record.getAttributeAsInt("ASM_End", -1)
                val asmStart = record.getAttributeAsInt("ASM_Start", -1)

                // length of ref segment and ASM segment
                // take absolute value of asmLen because if it is on the negative strand asmStart > asmEnd
                val refLen = record.end - record.start + 1
                val asmLen = abs(asmEnd - asmStart) + 1

                /* COVERAGE */

                // if we haven't encountered these chromosomes yet, add them to the map
                if (!baseCounterMap.containsKey(refChrom)) {
                    baseCounterMap[refChrom] = BaseCounter(TreeRangeSet.create<Int>(),
                        mutableListOf<Int>(), mutableListOf<Int>(), 0, 0, 0)
                }

                // add the ranges that are covered to the appropriate lists
                val refRange = getRange(record.start, record.end)

                if(baseCounterMap[refChrom]!!.rangesMapped.intersects(refRange)) {
                    myLogger.info("WARNING: Overlapping alignments found in chromosome $refChrom at position ${refRange.lowerEndpoint()}, ${refRange.upperEndpoint()}")
                }

                // add range to its respective contig set
                baseCounterMap[refChrom]!!.rangesMapped.add(refRange)

                /* INDEL AND SNPs */

                // we assume that these are haploid genotypes
                // can add support for diploids later

                val allele = record.genotypes[0].getAllele(0)
                val refAllele = record.reference

                // reference allele: add bases to numRefBases
                if (allele.isReference) {
                    baseCounterMap[refChrom]!!.ref += record.end - record.start + 1
                    // missing allele - represents N's
                } else if (allele.isNoCall) {
                    baseCounterMap[refChrom]!!.ns += record.end - record.start + 1
                    // not a symbolic allele: snp or indel
                } else if (!allele.isSymbolic) {
                    // SNP
                    if(refAllele.length() == 1 && allele.length() == 1) {
                        baseCounterMap[refChrom]!!.snps += 1
                        // indel
                    } else {
                        // positive is insertion, negative is deletion
                        // based on how BioKotlin's MafToVCF works, we shouldn't have equal lengths
                        val lengthDiff = allele.length() - refAllele.length()

                        if(lengthDiff > 0) {
                            baseCounterMap[refChrom]!!.ins.add(lengthDiff)
                        } else {
                            baseCounterMap[refChrom]!!.del.add(lengthDiff * -1)
                        }

                        // vcf indels require a padding base to be added to the beginning of the record
                        // this base could be either a SNP or reference
                        // so we want to check and add it to the appropriate count
                        if(refAllele.baseString[0] == allele.baseString[0]) {
                            baseCounterMap[refChrom]!!.ref += 1
                        } else {
                            baseCounterMap[refChrom]!!.snps += 1
                        }

                    }

                    // finally, for SNPs and indels add the number of Ns in the sample allele to the count
                    baseCounterMap[refChrom]!!.ns += allele.baseString.count{it == 'N'}

                }

            }

            // build GVCFStats objects
            val stats = baseCounterMap.map {
                GVCFStats(
                    sampleName,
                    it.key,
                    refChrLengths[it.key]!!,
                    it.value.snps,
                    it.value.ins.size,
                    it.value.del.size,
                    it.value.ns,
                    it.value.ins.sum(),
                    it.value.del.sum(),
                    it.value.ref.toDouble() / refChrLengths[it.key]!!,
                    it.value.rangesMapped.asRanges().map{
                            range -> range.upperEndpoint() - range.lowerEndpoint()}.sum().toDouble() / refChrLengths[it.key]!!,

                    if(it.value.ins.size > 0) {it.value.ins.sum().toDouble() / it.value.ins.size } else 0.0,
                    median(it.value.ins),
                    if(it.value.ins.size > 0) {it.value.ins.max() } else 0 ,

                    if(it.value.del.size > 0) { it.value.del.sum().toDouble() / it.value.del.size} else 0.0,
                    median(it.value.del),
                    if(it.value.del.size > 0) { it.value.del.max() } else 0
                )
            }.toMutableList()

            // add the aggregate data for the entire sample genome

            val totalRefLength = baseCounterMap.keys.sumOf { refChrLengths[it]!! }
            val totalIns = baseCounterMap.values.sumOf{it.ins.size}
            val totalDel = baseCounterMap.values.sumOf{it.del.size}

            stats.add(0, GVCFStats(
                sampleName,
                "ALL",
                totalRefLength,
                baseCounterMap.map{it.value.snps}.sum(),
                totalIns,
                totalDel,
                baseCounterMap.map{it.value.ns}.sum(),
                baseCounterMap.map{it.value.ins.sum()}.sum(),
                baseCounterMap.map{it.value.del.sum()}.sum(),
                baseCounterMap.map{it.value.ref}.sum().toDouble() / totalRefLength,
                baseCounterMap.map{it.value.rangesMapped.asRanges().map{
                        range -> range.upperEndpoint() - range.lowerEndpoint()}.sum()}.sum().toDouble() / totalRefLength,
                if(totalIns > 0) {baseCounterMap.map{it.value.ins.sum()}.sum().toDouble() / totalIns } else 0.0,
                median(baseCounterMap.flatMap{it.value.ins}),
                if(totalIns > 0) { baseCounterMap.flatMap{it.value.ins}.max() } else 0,
                if(totalDel > 0) {baseCounterMap.map{it.value.del.sum()}.sum().toDouble() / totalDel } else 0.0,
                median(baseCounterMap.flatMap{it.value.del}),
                if(totalDel > 0) { baseCounterMap.flatMap{it.value.del}.max() } else 0
            ))

            return stats

    }

    fun getHVCFStats(sample: String, rangeMap: Map<RefRangeInfo, HVCFStatsByRange>): List<HVCFStats> {
        val chromSet = rangeMap.keys.map { it.chrom }.toSet()

        val stats = chromSet.map { chrom ->
            val chromRanges = rangeMap.filterKeys { it.chrom == chrom }

            HVCFStats(
                sample,
                chrom,
                chromRanges.size,
                chromRanges.count { it.value.identicalToRef }
            )
        }.toMutableList()


        stats.add(0, HVCFStats(
            sample,
            "ALL",
            stats.map{it.refRangesWithHap}.sum(),
            stats.map{it.hapsIdenticalToRef}.sum()
        ))

        return stats

    }

    fun getHVCFStatsByRange(vcfFile: File): Pair<String, Map<RefRangeInfo, HVCFStatsByRange>> {

        VCFFileReader(vcfFile, false).use { reader ->
            val header = reader.fileHeader
            val records = reader.toList()

            return getHVCFStatsByRange(header, records)
        }
    }

    fun getHVCFStatsByRange(header: VCFHeader, records: List<VariantContext>): Pair<String, Map<RefRangeInfo, HVCFStatsByRange>> {

        val sample = header.genotypeSamples[0]

        val altHeaderLines = header.idHeaderLines.associate{ line ->
            try { Pair(line.id , line as VCFAltHeaderLine)
            } catch(exe: Exception) { Pair(line.id, null) }
        }.filterNot{ it.value == null }

        val refRangeMap = mutableMapOf<RefRangeInfo, HVCFStatsByRange>()
        val iterator = records.iterator()

        while(iterator.hasNext()) {
            val record = iterator.next()

            val haplotypeID = record.genotypes[0].getAllele(0).displayString.removeSurrounding("<", ">")

            if(!altHeaderLines.containsKey(haplotypeID)) {
                myLogger.warn("Haplotype $haplotypeID is not present in header.")
                continue
            }

            val headerLineFields = altHeaderLines[haplotypeID]!!.genericFields

            if(!headerLineFields.containsKey("Regions") || !headerLineFields.containsKey("RefRange")) {
                myLogger.warn("Header $haplotypeID is missing information.")
                continue
            }

            val ref = RefRangeInfo(headerLineFields["RefRange"]!!, record.contig, record.start, record.end)

            // from the regions, we can get the total length of the haplotype
            val sampleLength = headerLineFields["Regions"]!!.split(",").map{regionString ->
                val bounds = regionString.split(":")[1].split("-").map{it.toInt()}
                bounds[1] - bounds[0] + 1
            }.sum()

            val stats = HVCFStatsByRange(
                haplotypeID,
                sampleLength,
                haplotypeID == ref.id
            )

            refRangeMap[ref] = stats

        }

        return Pair(sample, refRangeMap)

    }


    /**
     * Standardizes the format of the ranges: closed on lower edge, open on top
     * This is mostly to make sure that concatenating two adjacent ranges
     * Produces one large range.
     * And flipping ranges where start > end
     */
    fun getRange(start: Int, end: Int): Range<Int> {
        return if(start > end) {
            Range.closed(end, start).canonical(DiscreteDomain.integers())
        } else {
            Range.closed(start, end).canonical(DiscreteDomain.integers())
        }
    }

    /**
     * Calculate the median value of a list of integers
     * Returns null if list is empty
     * **/
    fun median(values: List<Int>): Double {
        if (values.size == 0) { return 0.0 }

        val valuesSort = values.sorted()
        val mid = valuesSort.size / 2
        return if (valuesSort.size % 2 == 0) { (valuesSort[mid-1] + valuesSort[mid]).toDouble() / 2
        } else { valuesSort[mid].toDouble() }
    }
}