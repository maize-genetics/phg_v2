package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedReader
import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.vcf.VCFFileReader
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.parseALTHeader
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class provides metrics on the h.vcf file created from the imputation pipeline.
 * It takes an imputation h.vcf file as well as the h.vcf file from the sample that was
 * the parent of the short reads used in the imputation pipeline.  It also takes a bed file
 * that was used to create the ranges for the tiledb data storage.  It will output a file
 * that, for every range in the bed file, displays the haplotype IDs at that range in both
 * the imputation h.vcf as well as the parent sample h.vcf.  It will also display the read
 * counts for each haplotype ID to allow comparisons between the sample chosend at a given
 * range vs the read counts for the parent sample at that range.
 *
 * In addition, there are columns to note if the range was present in the imputation h.vcf
 * or in the parent sample h.vcf file, and a column that indicates if the 2 haplotypes (parent
 * and imputation) are the same.
 *
 * The output file will be tab-delimited and will be written to the output file specified
 * by the user.
 */

class ImputationMetrics : CliktCommand(help = "Impute best path(s) using read mappings.")   {
    private val myLogger = LogManager.getLogger(ImputationMetrics::class.java)
    val sampleHvcf by option(help = "Path to hvcf file loaded to tiledb for the Parent sample of your short reads. For example: If your reads are from P39, this would be the P39.h.vcf.")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }

    val imputationHvcf by option(help = "Path to the hvcf file created from the imputation pipeline. ")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }

    val bedFile by option(help = "Path to the bed file used to create ranges for the tiledb data storage. ")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }

    val chrom by option(help = "Comma separated list of chromosome to use for filtering output.  Only data for these chromosomes will be processed.")

    val outputFile by option(help = "Path to a file where the output data will be written. ")
        .required()

    val readMappingFiles by option(help = "Path to file containing a list of readMapping files created by PHG, 1 file per line.")
        .required()
        .validate { require(File(it).exists()) {"$it is not a valid file"} }


    override fun run() {
        processImputationResults(sampleHvcf, imputationHvcf, bedFile, chrom, outputFile, readMappingFiles)
    }

    fun processImputationResults(sampleHvcf: String, imputationHvcf: String, bedFile: String, chrom: String?, outputFile: String, readMappingFiles: String) {
        // create a List of chromosomes from the chrom string, which is a comma separated string
        val chromList = chrom?.split(",")?.toList()
        val parentChromPosToHapid = mutableMapOf<Position,String>()

        // Get just the name for the sampleHvcf file, and remove any h.vcf or hvcf extension
        val sampleHvcfName = File(sampleHvcf).name.replace(Regex("\\.hvcf$"), "").replace(Regex("\\.h.vcf$"), "")


        // Iterate through the sampleHvcf file and create a map of chrom/pos to hapid,
        // skipping positions on chromosomes not in the chromList
        // If the chromList is null or empty, process all chromosomes
        VCFFileReader(File(sampleHvcf),false).use { reader ->
            for (record in reader) {
                if (chromList == null || chromList.isEmpty() || chromList.contains(record.contig)) {
                    val alt = record.alternateAlleles[0].toString() // first value is the hapid
                    // from alt.string, strip off the leading and trailing "<>"
                    val hapid = alt.substring(1, alt.toString().length - 1)
                    parentChromPosToHapid[Position(record.contig, record.start)] = hapid
                }
            }
        }
        myLogger.info("Number of ranges in parentChromPosToHpaid: ${parentChromPosToHapid.size}")

        // Process the bedfile.  Only include ranges for chromosomes in the chromlist
        // If the chromslis is null or empty, include all ranges
        val bedRanges = mutableMapOf<Position, String>()
        bufferedReader(bedFile).useLines { lines ->
            lines.filter { line ->
                // If chromList is empty, allow all lines, otherwise filter based on chromList
                chromList.isNullOrEmpty() || chromList.any { chrom -> line.startsWith(chrom) }
            }.forEach { line -> // Process each filtered line
                val parts = line.split("\t")
                val start = parts[1].toInt() + 1 // BED is 0-based, VCF is 1-based
                val end = parts[2].toInt()
                val position = Position(parts[0], start) // Create Position object
                val refRange = "${parts[0]}:${start}-${end}" // Correctly formatted refRange string
                bedRanges[position] = refRange // Add entry to the map
            }
        }

        //Find which positions are in the bedfile, but not in the sampleHvcf (parentChromPosToHapid map)
        val sampleChromRanges = parentChromPosToHapid.keys.sorted()
        val sampleMissingRanges = bedRanges.keys.filter { !sampleChromRanges.contains(it) }

        // Read the imputation havcf file, get the altHeaders
        val reader2 = VCFFileReader(File(imputationHvcf),false)
        // THis is a map of Hapid -> AltHeaderMetaData
        // Can get the sampleName from the hapid, as well as the position.
        val altHeaderMap = parseALTHeader(reader2.header)
        reader2.close()
        myLogger.info("Num ranges in filtered bed: ${bedRanges.size}, num  ranges in sampleHvcf: ${sampleChromRanges.size}, num sampleHvcf missing ranges: ${sampleMissingRanges.size}")

        // Create a map of Position -> SampleName from the imputationHvcf file
        val wgsPosToHapid = mutableMapOf<Position, String>()
        val hapidToSampleName = mutableMapOf<String,String>()
        for (entry in altHeaderMap) {
            val hapid = entry.key
            val sampleName = entry.value.sampleName()
            // RefRange looks like this: RefRange=chr9:59390462-60165488, we get chr9:59390462-60165488 returned.
            // Create a Position object from the refRange, using the first part of the refRange as the contig up to the ":"
            // and from the ":" to the "-" as the start position.
            val refRange = entry.value.refRange //  Get the refRange from the AltHeaderMetaData
            val parts = refRange.split(":")
            val contig = parts[0]
            if (chromList.isNullOrEmpty() || chromList.contains(contig)) {
                val start = parts[1].split("-")[0].toInt()
                val position = Position(contig, start)
                wgsPosToHapid[position] = hapid
                hapidToSampleName[hapid] = sampleName
            }

        }
        val imputationMissingRanges = bedRanges.keys.filter { !wgsPosToHapid.keys.contains(it) }
        myLogger.info("Number of positions in filtered imputation file: ${wgsPosToHapid.keys.size}, number of Imputation missing ranges: ${imputationMissingRanges.size}")


        // Create map of haplotype Id to read counts
        val hapidCounts = readCountsFromFiles(readMappingFiles)

        // Create header line for tab-delimited output file
        val headerLine = "refRange\tHapotypeId\tImputationSampleName\tIn_ImputationHvcf\tIn_${sampleHvcfName}\t${sampleHvcfName}_HaplotypeId\tImputationSampleRdCount\t${sampleHvcfName}SampleRdCt\tHapidsMatch\n"

        bufferedWriter(outputFile).use {writer ->
            writer.write(headerLine)
            // Iterate through the bedChr9Ranges, and for each position, get the hapid from the wgsPosToHapid map
            // and the sampleName from the hapidToSampleName map.  Get the readCount for the hapid from the hapidCounts map.
            // Get the readCount for the P39 haplotype from the hapidCounts map.
            // Write the refRange, sampleName, inImputationHvcf, inP39Hvcf, P39RdCount, SampleRdCount to the file.
            // Zack wants a histogram of the difference between SampleRdCount and P39RdCount
            // create a list, then will write it to a file below

            for (pos in bedRanges.keys.sorted()) {
                val refRange = bedRanges[pos]
                val hapid = wgsPosToHapid.getOrDefault(pos, "NA")
                val inImpuptation = if (hapid == "NA") "No" else "Yes"
                val inP39 = if (parentChromPosToHapid.containsKey(pos)) "Yes" else "No"
                // hapid always has a value: NA if not present, and NA won't be in hapidToSampleName, so this works
                val sampleName = hapidToSampleName.getOrDefault(hapid, "NA")
                // get readCount for the hapid
                val sampleRdCount = hapidCounts.getOrDefault(hapid, 0)
                // get readCount for parentSampleHapid
                val parentSampleHapid = parentChromPosToHapid.getOrDefault(pos, "NA")
                val parentSampleRdCt = if (parentSampleHapid == "NA") 0 else hapidCounts.getOrDefault(parentSampleHapid, 0)
                // check if the hapid and parentSampleHapid are the same
                val hapidMatch = if (hapid != "NA" && parentSampleHapid != "NA" && hapid == parentSampleHapid) "true" else "false"

                writer.write("$refRange\t$hapid\t$sampleName\t$inImpuptation\t$inP39\t$parentSampleHapid\t$sampleRdCount\t$parentSampleRdCt\t$hapidMatch\n")
            }
        }
    }

    // function to read the readMapping files and create a map of Hapid -> read counts
    fun readCountsFromFiles(readMappingFiles:String): Map<String, Int> {
        // The readMappingFiles parameter is file containing full path to a read mapping file, one per line.
        // Each individual file in this list has headers, then is a tab-delimited file with the following columns:
        // HapIds (comma-separated list)
        // count
        // Read each file, and for each file:
        //  - split up the hapids, and create a map of haplotypeID->readCount
        //  - For each hapid in the split list, add the readCount to the map for that hapid
        // Create a map of haplotypeID->readCount
        val hapidCounts = mutableMapOf<String, Int>()
        bufferedReader(readMappingFiles).useLines { lines ->
            lines.forEach { line ->
                // Skip blank lines
                if (line.isBlank()) return@forEach
                myLogger.info("Processing file: $line")
                bufferedReader(line).useLines { fileLines ->
                    fileLines.forEach { fileLine ->
                        // If line begins with # or HapIds then skip
                        if (!fileLine.startsWith("#") && !fileLine.startsWith("HapIds")) {
                            val parts = fileLine.split("\t")
                            val hapids = parts[0].split(",")
                            val count = parts[1].toInt()
                            hapids.forEach { hapid ->
                                hapidCounts[hapid] = hapidCounts.getOrDefault(hapid, 0) + count
                            }
                        }
                    }
                }
            }
        }
        return hapidCounts
    }


}