package net.maizegenetics.phgv2.pathing.ropebwt

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.mutuallyExclusiveOptions
import com.github.ajalt.clikt.parameters.groups.required
import com.github.ajalt.clikt.parameters.groups.single
import com.github.ajalt.clikt.parameters.options.convert
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.cli.logCommand
import net.maizegenetics.phgv2.pathing.AlignmentUtils
import net.maizegenetics.phgv2.pathing.KeyFileData
import net.maizegenetics.phgv2.pathing.ReadInputFile
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.BufferedReader
import java.io.InputStreamReader

/**
 * MapReadsVcf will map VCF files to a pangenome indexed by ropeBWT3
 * This will create a standard read mapping file that can be used for path finding
 * Note that this will take some time to run.  Our internal tests show that WGS files can take about 20 minutes to run.
 * The VCF file is first used to create a fastq file by extracting a 200 base segment from the Reference FASTA.
 */
class MapReadsVcf : CliktCommand(help="Map VCF to a pangenome using ropeBWT3") {

    private val myLogger = LogManager.getLogger(MapReadsVcf::class.java)

    val index by option(help = "The full path of the ropebwt3 index file.")
        .required()

    val readInputFiles: ReadInputFile by mutuallyExclusiveOptions<ReadInputFile>(
        option("--key-file", help = "Name of tab-delimited key file.  Columns for samplename and filename are required.  If using paired end fastqs, a filename2 column can be included. A value must be entered for either --key-file or --read-files.").convert{ ReadInputFile.KeyFile(it) },
        option("--read-files", help = "Comma separated list of fastq files for a single sample.  Either 1(for single end) or 2(for paired end) files can be input at a time this way.  Any more and an error will be thrown.").convert{ ReadInputFile.ReadFiles(it) }
    ).single().required()

    val outputDir by option("-o", "--output-dir", help = "Name for output ReadMapping file Directory (Required)")
        .required()

    val hvcfDir by option(help = "Directory for the haplotype VCF files. These hvcfs will be used to filter the hits to only those that hit a single referenceRange.")
        .required()

    val threads by option(help = "Number of threads to use.")
        .int()
        .default(5)

    val minMemLength by option(help = "Minimum length of a match to be considered a match.")
        .int()
        .default(148)

    val maxNumHits by option(help = "Number of hits to report.  Note ropebwt can hit more than --max-num-hits but any alignment hitting more haplotypes than this will be ignored.")
        .int()
        .default(50)

    val condaEnvPrefix by option (help = "Prefix for the conda environment to use.  If provided, this should be the full path to the conda environment.")
        .default("")

    val maxStart by option(help = "Maximum start position for a read to be considered a match. Any alignments with a start above this will be ignored.")
        .int()
        .default(0)

    val minEnd by option(help = "Minimum end position for a read to be considered a match. Any alignments with an end below this will be ignored.")
        .int()
        .default(70)

    data class VcfRecord(
        val chrom: String,
        val pos: Int,
        val id: String,
        val ref: String,
        val alt: String,
        val genotypes: List<String>,
        var fastq: String
    )

    data class VcfResult(
        val header: List<String>,
        val records: List<VcfRecord>
    )

    override fun run() {
        logCommand(this)

        myLogger.info("Building the Graph")

        val graph = HaplotypeGraph(hvcfDir)

        val hapIdToRefRangeMap = graph.hapIdToRefRangeMap()

        myLogger.info("Mapping reads to pangenome")

        mapAllReadFiles(index, readInputFiles.getReadFiles(), outputDir, threads, minMemLength, maxNumHits, condaEnvPrefix, hapIdToRefRangeMap, maxStart, minEnd)
    }

    /**
     * Function to map all the VCF file in the keyFileDataEntries to the index
     */
    fun mapAllReadFiles(index: String, keyFileDataEntries: List<KeyFileData>, outputDir: String, threads: Int,
                        minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String,
                        hapIdToRefRangeMap : Map<String, List<ReferenceRange>>, maxStart: Int, minEnd: Int) {
        //Loop through the keyFileDataEntries and map the reads
        //file1 is fasta
        //file2 is VCF
        val readNameToFileMap = mutableMapOf<String,MutableList<String>>()
        for(readFile in keyFileDataEntries) {
            val fileList = readNameToFileMap[readFile.sampleName]?: mutableListOf()
            mapFastaReadFile(index, readFile.sampleName, readFile.file1, readFile.file2, outputDir, threads, minMemLength, maxNumHits, condaEnvPrefix, hapIdToRefRangeMap, maxStart, minEnd)
            //fileList.add(outputFile1)
            readNameToFileMap[readFile.sampleName] = fileList
        }
        exportPathKeyFile(outputDir, readNameToFileMap)
    }

    /**
     * Read in FASTA file that has been generated from samtools faidx
     */
    fun parseFastaOutput(fastaFilePath: String): Map<String, String> {
        val sequenceMap = mutableMapOf<String, String>()
        val lines = File(fastaFilePath).readLines()
        var currentId = ""
        val seqBuilder = StringBuilder()

        for (line in lines) {
            if (line.startsWith(">")) {
                // Save previous
                if (currentId.isNotEmpty()) {
                    sequenceMap[currentId] = seqBuilder.toString()
                    seqBuilder.clear()
                }
                currentId = line.removePrefix(">").trim()
                if (sequenceMap.containsKey(currentId)) {
                    println("Warning: Duplicate header detected: $currentId")
                }
            } else {
                seqBuilder.append(line.trim())
            }
        }
        if (currentId.isNotEmpty()) {
            sequenceMap[currentId] = seqBuilder.toString()
            //println("$currentId ${sequenceMap[currentId]}")
        }
        return sequenceMap
    }

    /**
     * add fasta to VCF records
     */
    fun mapFastaToVcfRecords(
        vcfRecords: List<VcfRecord>,
        fastaSequences: Map<String, String>,
        flankSize: Int = 100
    ): List<VcfRecord> {
        return vcfRecords.map { record ->
            val start = record.pos - flankSize
            val end = record.pos + flankSize
            val key = "${record.chrom}:${start}-${end}"

            val seq = fastaSequences[key] ?: ""
            record.copy(fastq = seq)
        }
    }

    /**
     * Read the VCF file and extract the header, ref, alt, position, and genotypes
     */
    fun readVcfAndExtractRegions(filePath1:String, filePath2: String, regionsPath: String, fqPath: String): VcfResult {
        val file = File(filePath2)
        var header: List<String> = emptyList()
        val records = mutableListOf<VcfRecord>()
        val regionKeys = mutableListOf<String>()
        var count = 0

        if (!file.exists()) {
            println("File not found: $filePath2")
            return VcfResult(header, records)
        }

        file.forEachLine { line ->
            when {
                line.startsWith("##") -> {
                    // Meta-information header lines
                    // println("Meta: $line")
                }
                line.startsWith("#") -> {
                    // Column header line
                    header = line.split("\t").drop(9)
                    myLogger.info("header " + header)
                }
                else -> {
                    // Variant data lines
                    val fields = line.split("\t")
                    if (fields.size >= 8) {
                        // only use gt entries
                        // convert 0/1 to allele
                        // I don't see support for missing or ambiguous codes so use reference
                        val chrom = fields[0]
                        val pos = fields[1].toInt()
                        val id = fields[2]
                        val ref = fields[3]
                        val alt = fields[4]
                        count++
                        val genotypes = fields.drop(9).map { call ->
                            val gt = call.split(":").getOrElse(0) { "." }
                            when (gt) {
                                "0/0", "./.", ".", "0/1" -> ref
                                "1/1" -> alt
                                else -> {
                                    print("undefined gt $gt")
                                    ref
                                }
                            }
                        }
                        val region = "$chrom:${pos - 100}-${pos + 100}"
                        regionKeys.add(region)
                        records.add(
                            VcfRecord(chrom, pos, id, ref, alt, genotypes, fastq = "")
                        )
                    }
                }
            }
        }
        myLogger.info("readVcf found $count entries")
        // Write regions to a temp file
        //val regionsFile = File("/project/genolabswheatphg/clay_phg/output/mapping/fastq/regions.txt")
        val regionsFile = File(regionsPath)
        if (regionsFile.exists()) {
            regionsFile.delete()
        }
        regionsFile.printWriter().use { out ->
            regionKeys.forEach { out.println(it) }
        }

        // Fetch sequences all at once
        myLogger.info("regionsFile: ${regionsFile.absolutePath}")
        val command = listOf("samtools", "faidx", "-r", regionsFile.absolutePath, filePath1, "-o", fqPath)
        val process = ProcessBuilder(command)
            .redirectErrorStream(true)
            .start()
        val exitCode = process.waitFor()
        if (exitCode != 0) {
            val err = process.inputStream.bufferedReader().readText()
            throw RuntimeException("samtools failed: $err")
        }

        myLogger.info("Processed ${records.size} variants with batched samtools call.")
        return VcfResult(header, records)
    }

    /**
     * Function to map a VCF generated read file to the index and write the read mapping to the outputFile
     */
    fun mapFastaReadFile(index: String, sampleName: String, readFile1: String, readFile2: String, path: String, threads: Int,
                          minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String,
                          hapIdToRefRangeMap: Map<String,List<ReferenceRange>>, maxStart: Int, minEnd: Int) {
        myLogger.info("Mapping reads in VCF $readFile2 using FASTA $readFile1 to $index")
        val regionsPath = "/project/genolabswheatphg/clay_phg/output/mapping/fastq/regions.txt"
        val fqPath = "/project/genolabswheatphg/clay_phg/output/mapping/fastq/regions.fq"
        val results = readVcfAndExtractRegions( filePath1= readFile1, filePath2 = readFile2, regionsPath, fqPath)
        val readFile4 = "$path/fastq/regions.fq"
        val fastaMap = parseFastaOutput(readFile4)
        val updatedRecords = mapFastaToVcfRecords(results.records, fastaMap)
        myLogger.info("Done adding fasta to VCF records")
        val qual = CharArray(201) { 'I' }.concatToString()
        val dir = File(path)
        if(!dir.exists()) {
            dir.mkdirs()
        }

        results.header.forEachIndexed {index1, sample ->
            myLogger.info("processing sample $index1 $sample")
            val sb = StringBuilder()
            updatedRecords.forEach { record ->
                var id = record.id
                val ref = record.ref
                val alt = record.alt
                val genotypes = record.genotypes
                val cleanSequence = record.fastq
                if (!id.contains(Regex("[a-zA-Z]"))) {
                    id = record.chrom + "_" + record.pos
                }
                try {
                    if (ref[0] == cleanSequence[100]) {
                        //println("Good: $id $ref")
                    } else {
                        println("Bad: $id $cleanSequence ${ref[0]} ${alt[0]} ${cleanSequence[100]}")
                    }
                    val gt = genotypes[index1]
                    val newStr = cleanSequence.toCharArray().also { it[100] = gt[0] }.concatToString()
                    sb.appendLine("@$id\n$newStr\n+\n$qual\n")
                } catch (e: Exception) {
                    println("Error: ${e.message}")
                }
            }
            val safeSample = sample.replace("/","_")
            val file = File("$path/fastq/" + safeSample + ".fq")
            if (file.exists()) {
                file.delete()
            }
            file.writeText(sb.toString())
            myLogger.info("Created file: $sample")
        }
        myLogger.info("Done creating fastq files from VCF")

        //now read from created fastq file
        results.header.forEach { filename ->
            val safeFilename = filename.replace("/","_")
            val readFile3 = path + "/fastq/" + safeFilename + ".fq"
            myLogger.info("Setup mapping to $readFile3")
            val bedFileReader =
                    setupMappingProcess(index, readFile3, threads, minMemLength, maxNumHits, condaEnvPrefix)
            val readMapping =
                    createReadMappingsForFileReader(bedFileReader, maxNumHits, hapIdToRefRangeMap, maxStart, minEnd)
            bedFileReader.close()
            //val outputFile2 = outputFile + "/" + filename + "_readMapping.txt"
            val outputFile2 = "$path/" + safeFilename + "_readMapping.txt"
            myLogger.info("Writing read mapping to $outputFile2")
            AlignmentUtils.exportReadMapping(outputFile2, readMapping, sampleName, Pair(readFile1, ""))
        }
    }

    /**
     * Function to setup the ropebwt3 mem process and pass a BufferedReader for use by the rest of the program
     * //time ../ropebwt3/ropebwt3 mem -t40 -l148 -p50 /workdir/zrm22/phgv2/ropeBWT/fullASMTests/phg_ASMs.fmd /workdir/zrm22/phgv2/ropeBWT/Reads/B97_HVMFTCCXX_L7_1.clean.fq.gz > B97_1_fullASM_pos_matches2NM.bed
     *
     */
    fun setupMappingProcess(index: String, readFile: String, threads: Int, minMemLength: Int, maxNumHits: Int, condaEnvPrefix: String): BufferedReader {
        val prefixArg = if(condaEnvPrefix.isNotBlank()) {
            Pair("-p",condaEnvPrefix)
        }
        else {
            Pair("-n", "phgv2-conda")
        }

        val ropebwt3Process = ProcessBuilder("conda","run",prefixArg.first,prefixArg.second,"ropebwt3", "mem", "-t$threads", "-l$minMemLength", "-p$maxNumHits", index, readFile)
            .redirectError(ProcessBuilder.Redirect.INHERIT)
            .start()

        return BufferedReader(InputStreamReader(ropebwt3Process.inputStream))
    }

    /**
     * Function to create a readMapping Map from a BufferedReader of a ropebwt3 mem output
     * This will work directly from a file but this command is not setup to do that yet.
     */
    fun createReadMappingsForFileReader(
        bedFileReader: BufferedReader,
        maxNumHits: Int,
        hapIdToRefRangeMap: Map<String, List<ReferenceRange>>,
        maxStart: Int,
        minEnd: Int
    ): MutableMap<List<String>, Int> {
        var currentLine = bedFileReader.readLine()
        val tempMems = mutableListOf<MEM>()
        val readMapping = mutableMapOf<List<String>, Int>()
        while (currentLine != null) {
            if(currentLine.isEmpty()) {
                currentLine = bedFileReader.readLine()
                continue
            }
            val alignmentParsed = RopeBWTUtils.parseStringIntoMem(currentLine)
            if (tempMems.isNotEmpty() && tempMems[0].readName != alignmentParsed.readName) {
                //write out the tempMems
                processMemsForRead(tempMems, readMapping, maxNumHits,hapIdToRefRangeMap, maxStart, minEnd)
                tempMems.clear()
            }
            tempMems.add(alignmentParsed)
            currentLine = bedFileReader.readLine()
        }

        if(tempMems.isNotEmpty()) {
            processMemsForRead(tempMems, readMapping, maxNumHits, hapIdToRefRangeMap, maxStart, minEnd)
        }

        return readMapping
    }

    /**
     * Function to process the mems for a single read and add them to the readMapping
     * This will filter things based on the maxStart and minEnd positions, then  maxNumHits and retain the mems that are longest
     * Because MEMs are Maximal Exact Matches, the longest MEMs are the best hits
     */
    fun processMemsForRead(tempMems: List<MEM>, readMapping: MutableMap<List<String>, Int>, maxNumHits: Int,
                           hapIdToRefRangeMap: Map<String, List<ReferenceRange>>, maxStart: Int, minEnd: Int) {
        val posFilteredMems = tempMems.filter { it.readStart <= maxStart && it.readEnd >= minEnd }
        if(posFilteredMems.isEmpty()) {
            return
        }
        //get the longest hits
        val maxLength = posFilteredMems.maxOf { it.readEnd - it.readStart }
        //remove any hits that are not the longest
        val bestHits = posFilteredMems.filter { it.readEnd - it.readStart == maxLength }

        val totalNumHits = bestHits.sumOf { it.numHits }

        if(totalNumHits <= maxNumHits) {
            //if the total number of hits is less than the maxNumHits, Filter the haps down to a single ref range and then add all those to the readMapping
            //First turn the hapIds into a set then back to a list and sort so it will be consistent
            val bestHapIds = bestHits.flatMap { it.listMemHits.map{hits -> hits.contig} }.toSet()

            val filteredBestHapIds = filterToOneReferenceRange(bestHapIds, hapIdToRefRangeMap)

            val hapIdsHit = filteredBestHapIds.sorted()
            readMapping[hapIdsHit] = (readMapping[hapIdsHit]?:0) + 1
        }
    }

    /**
     * Function to filter the hapIds to only those that hit a single reference range
     * This will pick the highest occurring reference range and if there is a tie it will choose the first one.
     */
    fun filterToOneReferenceRange(hapIds: Set<String>, hapIdToRefRangeMap: Map<String, List<ReferenceRange>>) : List<String> {
        //figure out which reference range is hit the most
        val refRangeCounts = hapIds.map { hapIdToRefRangeMap[it] }
            .filterNotNull()
            .flatten()
            .groupingBy { it }
            .eachCount()

        val maxRefRange = refRangeCounts.maxByOrNull { it.value }?.key
        //Then filter out any hapIds that don't have that reference range
        return hapIds.filter { hapIdToRefRangeMap[it]?.contains(maxRefRange)?:false }
    }

    /**
     * Function to export the path key file for the readMapping files that have been processed.
     */
    fun exportPathKeyFile(outputDir: String, readNameToFileMap: Map<String, List<String>>) {
        val pathKeyFile = "$outputDir/pathKeyFile.txt"
        myLogger.info("Writing pathKeyFile to $pathKeyFile")
        bufferedWriter(pathKeyFile).use {writer ->
            writer.write("sampleName\tfilename\n")
            for ((sampleName, fileList) in readNameToFileMap) {
                for (file in fileList) {
                    writer.write("$sampleName\t$file\n")
                }
            }
            writer.close()
        }
    }
}
