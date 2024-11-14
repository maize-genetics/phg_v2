package net.maizegenetics.phgv2.pathing

import biokotlin.util.bufferedWriter
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.flag
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.int
import com.github.ajalt.clikt.parameters.types.long
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap
import it.unimi.dsi.fastutil.longs.LongOpenHashSet
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.pathing.AlignmentUtils.Companion.buildHaplotypeGraph
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.time.DurationUnit
import kotlin.time.measureTimedValue

data class IndexedKmerData(val chr: String, val pos: UInt, val gameteLong0 : ULong, val gameteLong1 : ULong)
class BuildRamEfficientKmerIndex: CliktCommand(help="Create a kmer index for a HaplotypeGraph. By default the file will be written " +
        "to <hvcfDir>/kmerIndex.txt") {

    private val myLogger = LogManager.getLogger(BuildKmerIndex::class.java)

    val dbPath by option(help = "Tile DB URI")
        .required() //Needs to be required now due to the agc archive

    val indexFile by option(help = "The full path of the kmer index file. Default = <hvcf-dir>/kmerIndex.txt")
        .default("")

    val discardFile by option(help = "The full path of the discard set file. If Left blank no file will be written out.")
        .default("")

    val maxHaplotypeProportion by option("-p", "--maxHapProportion", help = "only kmers mapping to less than or " +
            "equal to maxHapProportion of haplotypes in a reference range will be retained.")
        .double()
        .default(0.75)

    val hashMask by option("-m", "--hashMask", help = "with hashFilter, used to mask kmers for filtering. " +
            "Default uses only the last kmer nucleotide. Only change this if you know what you are doing.")
        .long()
        .default(3)

    val hashFilterValue by option("-f", "--hashFilter", help = "Only hashes that pass the filter" +
            " ((hashValue and hashMask) == hashFilter) will be considered. Do not change this value unless you know " +
            "what you are doing.")
        .long()
        .default(1)

    val hvcfDir by option("--hvcf-dir", help = "Path to directory holding hVCF files. Data will be pulled directly from these files instead of querying TileDB")
        .required()//Todo: make this optional by adding .default("")

    val maxArgLength by option(help="The maximum argument length for a call to agc." +
            "If you get an error caused by a call to agc being too long try reducing this value.")
        .int()
        .default(200000)

    val noDiagnostics by option("-n", "--no-diagnostics", help = "Flag that will suppress writing of diagnostics.").flag()

    val initialKeepSize by option(help="Initial size of the keep map.  Default is 500,000,000")
        .int()
        .default(500_000_000)

    private val refrangeToAdjacentHashCount = mutableMapOf<ReferenceRange, Int>()
    private var runDiagnostics = true
    override fun run() {

        //build the haplotypeGraph
        myLogger.info("Start of BuildRamEfficientKmerIndex...")
        val graph = buildHaplotypeGraph(hvcfDir)

        val (hashToHapidMap, discardSet) = processGraphKmers(graph, dbPath, maxHaplotypeProportion,  hashMask, hashFilterValue, initialKeepSize)

        val kmerIndexFilename = if (indexFile == "") "${hvcfDir}/kmerIndex.txt" else indexFile


        //save the kmerIndex
        saveKmerHashesAndHapids(graph, kmerIndexFilename, hashToHapidMap)

        if(discardFile.isNotEmpty()) {
            bufferedWriter(discardFile).use { writer ->
                for (element in discardSet) {
                    writer.write("$element\n")
                }
            }
        }

        if (noDiagnostics) {
            runDiagnostics = false
            myLogger.info("BuildKmerIndex: Diagnostic output will not be written because the --no-diagnostic flag was set.")
        }
        else {
            writeDiagnostics(refrangeToAdjacentHashCount)
        }

    }


    // This is derived from the function processGraphKmers in BuildKmerIndex.kt
    // This implementation uses a mutableMap, which should be able to
    // grow as needed.  The bitmaps that contain the gamete information will take up
    // less room and need less storage.
    fun processGraphKmers(graph: HaplotypeGraph, dbPath: String, maxHaplotypeProportion: Double=.75,
                          hashMask: Long = 3, hashFilterValue:Long = 1, initialKeepSize: Int = 500_000_000) : Pair<Map<Long,IndexedKmerData>, LongOpenHashSet> {

        // keepMap is a map of hapidHash -> IndexedKmerData
        val keepMap = mutableMapOf<Long,IndexedKmerData>()
        //discardSet is a Set of hashes
        val discardSet = LongOpenHashSet(initialKeepSize)
        val startTime = System.nanoTime()

        // This returns a SortedSEt of <SampleGamete> - should it be turned into an ordered
        // list of just sample names?  This list is needed to index into the bit arrays
        val sampleGametes = graph.sampleGametesInGraph()
        val samples = sampleGametes.map { it.name }.sorted()

        // TODO - implement the rest of this!

        return Pair(keepMap, discardSet)

    }

    private fun saveKmerHashesAndHapids(
        graph: HaplotypeGraph,
        kmerIndexFilename: String,
        hashToHapidMap: Map<Long, IndexedKmerData>
    ) {
        TODO("Not yet implemented")
    }

    private fun writeDiagnostics(refrangeToAdjacentHashCount: MutableMap<ReferenceRange, Int>) {
        TODO("Not yet implemented")

    }

}