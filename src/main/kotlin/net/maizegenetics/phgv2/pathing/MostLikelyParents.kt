package net.maizegenetics.phgv2.pathing

import com.google.common.collect.Multimap
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.api.ReferenceRange
import net.maizegenetics.phgv2.api.SampleGamete

/**
 * Given a set of read mappings, finds the most likely parents in a stepwise fashion. The parent with the most reads
 * is selected as a likely parent. Then the the subset of reads not mapping to that parent is created.
 * The parent with the most reads from that subset is then added to the list of likely parents.
 * The process is repeated until the number of likely parents = maxParents, coverage is equal to or greater than minCoverage,
 * or until all parents have been added to the list, whichever occurs first. Coverage is calculated as the number
 * of reads mapping to any of the likely parents divided by the total number of reads.
 */
class MostLikelyParents
/**
 * A constructor that takes a [HaplotypeGraph].
 */(hapGraph: HaplotypeGraph) {
    //list of parents, which are the samples (eventually sampleGametes) present in the HaplotypeGraph
    val myParentList: List<SampleGamete>
    //a map of ReferenceRange -> map (parent -> hapid)
    val myParentToHapidMapByRefRange: Map<ReferenceRange, Map<SampleGamete,String>>

    init {
        myParentList = getSampleGametes(hapGraph)
        myParentToHapidMapByRefRange = refRangeToMapOfSampleToHapid(hapGraph)
    }

    /**
     * A constructor that takes a list of parents and for each ReferenceRange, a map of parent name to its hapid.
     * This constructor is intended mainly for unit testing.
     */
//    constructor(parentList: List<String>, parentHapidMap: Map<ReferenceRange, Map<String,String>>) {
//        myParentList = parentList
//        myParentToHapidMapByRefRange = parentHapidMap
//    }

    /**
     * Finds the most likely parents for a set of hapid counts from read mappings.
     */
    fun findMostLikelyParents(refRangeToHapIdSetCounts: Multimap<ReferenceRange, HapIdSetCount>, maxParents: Int, minCoverage: Double) : List<Pair<SampleGamete, Int>> {
        //convert refRangeToHapIdSetCounts to a Map<ReferenceRange,List<HapIdSetCount>>
        //use only the counts for ranges that are present in myParentToHapidMap
        var filteredCounts = refRangeToHapIdSetCounts.entries().filter { (refrange, _) -> myParentToHapidMapByRefRange.keys.contains(refrange) }
            .groupBy({it.key},{it.value})
        val bestParentList = mutableListOf<Pair<SampleGamete, Int>>()
        var iteration = 0
        var coverage = 0.0
        val totalCount = filteredCounts.map {(_,setCounts) -> setCounts.sumOf { it.count } }.sum()
        var cumulativeCount = 0

        while (iteration < maxParents && coverage < minCoverage) {
            iteration++
            var bestParent = SampleGamete("none", 0)
            var highestCount = 0
            for (parent in myParentList) {
                //get count for this parent
                val parentCount = //for each refrange count parent use
                    filteredCounts.keys.sumOf { refrange ->
                        //for each refrange count parent use
                        val hapid = myParentToHapidMapByRefRange[refrange]!![parent]
                        if (hapid == null) 0 else
                            filteredCounts[refrange]!!.filter { it.hapIdSet.contains(hapid) }.sumOf { it.count }
                    }
                if (parentCount > highestCount) {
                    highestCount = parentCount
                    bestParent = parent
                }
            }

            bestParentList.add(Pair(bestParent, highestCount))
            cumulativeCount += highestCount
            coverage = cumulativeCount / totalCount.toDouble()

            //rebuild the filteredList using only hapIdSetCounts that do not contain the best parent
            //this will be used for the next round
            filteredCounts = filteredCounts.keys.map { refrange ->
                //for each refrange count parent use
                val hapid = myParentToHapidMapByRefRange[refrange]!![bestParent]
                val filteredList = filteredCounts[refrange]!!.filter { !it.hapIdSet.contains(hapid) }
                Pair(refrange, filteredList)
            }.filter { it.second.isNotEmpty() }.toMap()
        }

        return bestParentList
    }


    companion object {

        data class HapIdSetCount(val hapIdSet : Set<String>, val count : Int)
        fun getSampleGametes(hapGraph: HaplotypeGraph) : List<SampleGamete> {
            val sampleSet = mutableSetOf<SampleGamete>()
            for (refrange in hapGraph.ranges()) {
                hapGraph.hapIdToSampleGametes(refrange).values.forEach { sampleSet.addAll(it) }
            }
            return sampleSet.sorted()
        }

        fun refRangeToMapOfSampleToHapid(hapGraph: HaplotypeGraph) : Map<ReferenceRange, Map<SampleGamete,String>> {
            val parentMap = mutableMapOf<ReferenceRange, Map<SampleGamete, String>>()
            for (refrange in hapGraph.ranges()) {
                val sampleSet = hapidToSampleGametes(hapGraph, refrange).map { it.value }.flatten().toSet()
                val rangeParentMap = sampleSet.associateWith { hapGraph.sampleToHapId(refrange, it) }
                parentMap.put(refrange, rangeParentMap)
            }
            return parentMap
        }

        fun refRangeToMapOfSampleToHapid(hapGraph: HaplotypeGraph, sampleSet: List<SampleGamete>) : Map<ReferenceRange, Map<SampleGamete,String>> {
            val parentMap = mutableMapOf<ReferenceRange, Map<SampleGamete, String>>()
            for (refrange in hapGraph.ranges()) {
                //TODO deal with case where sample is not in this reference range
                val rangeParentMap = sampleSet.associateWith { hapGraph.sampleToHapId(refrange, it) }
                parentMap.put(refrange, rangeParentMap)
            }
            return parentMap
        }

        fun hapidToSampleGametes(graph: HaplotypeGraph, refrange: ReferenceRange): Map<String, List<SampleGamete>> {
            return graph.hapIdToSamples(refrange).entries.map { (hapid, nameList) ->
                hapid to nameList.map{SampleGamete(it)}
            }.toMap()
        }
    }
}

/**
 * This plugin finds the most likely parents for each of a series of read mapppings specified by either a read method
 * or a key file. The parents are selected in a stepwise manner by finding the parent with the most mapped reads not
 * mapped to previously selected parents. Additional parents are selected until maxParents have been selected,
 * all available parents are selected, or minCoverage is reached. minCoverage is calculated as the number of reads
 * mapping to any of the likely parents divided by the total number of reads.
 *
 * One [HaplotypeGraph] is required as input.
 *
 * The plugin gets DB connection parameters from the ParameterCache, so expects that to be populated.
 *
 * The output file has headers: sample_name, total_reads, parent, order, number_of_reads, cumulative_reads
 */

//TODO convert to Clickt class
/*
class LikelyParentsPlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    private val myLogger = Logger.getLogger(LikelyParentsPlugin::class.java)

    private var outputFile = PluginParameter.Builder("outFile", null, String::class.java)
        .description("The file to which the output will be written. If the file exists, an error will be thrown.")
        .outFile()
        .required(true)
        .build()

    private var overwriteOutput = PluginParameter.Builder("overwrite", false, Boolean::class.javaObjectType)
        .description("If this is set to true and the output file exists, it will be overwritten without warning.")
        .build()

    private var keyFile = PluginParameter.Builder("keyFile",null, String::class.java)
        .description("KeyFile file name.  Must be a tab separated file using the following headers:\n" +
                "SampleName\tReadMappingIds\n" +
                "Additional headers will be ignored. ReadMappingIds need to be comma separated for multiple values. \n" +
                "If a keyFile is not supplied then all of the read mappings for the supplied method will be run.")
        .required(false)
        .inFile()
        .build()

    private var readMethod = PluginParameter.Builder("readMethod", null, String::class.java)
        .description("A read method. Either a read method or a keyfile containing mapping ids must be supplied.\n" +
                "If both are supplied, the keyFile will be used and the readMethod will be ignored.")
        .build()

    private var maxNumberOfParents = PluginParameter.Builder("maxParents", Int.MAX_VALUE, Int::class.javaObjectType)
        .description("The maximum number of parents to be selected.")
        .build()

    private var minCoverage = PluginParameter.Builder("minCoverage", 1.00, Double::class.javaObjectType)
        .description("No more parents will be selected when the proportion of total reads mapped to the " +
                "selected parents is greater than or equal to this number.")
        .build()

    private lateinit var readMappingDecoder: ReadMappingDecoder

    override fun pluginDescription(): String {
        return "The LikelyParentsPlugin finds the parents with the highest marginal read count fit in a forward stepwise " +
                "fashion for each sample or keyfile entry. It starts by finding the parent with the most reads. " +
                "Then, it finds the parent with the most reads in the set of reads not that did not map to the first parent. " +
                "It proceeds to find other parents in a similar stepwise manner, stopping when the first " +
                "limiting condition is reached. Possible limiting conditions are maximum number of parents " +
                "and proportion of total reads fit. Each parent selected together with the marginal number " +
                "of reads mapped to it are output. The plugin requires that the pipeline argument -configParameters " +
                "be used for the database config file."
    }

    override fun processData(input: DataSet?): DataSet? {
        //check argument requirements
        require(keyFile() != null || readMethod() != null) {"Both keyFile and readMethod were null. A value must be supplied for one of them."}
        require(input != null) {"A HaplotypeGraph must be supplied as input"}
        if (!overwriteOutput()) require(!File(outputFile()).exists()) {"The file ${outputFile()} already exists and will not be overwritten."}
        val dataList = input.getDataOfType(HaplotypeGraph::class.java)
        require(dataList.size == 1) {"Exactly one HaplotypeGraph must be supplied as input"}

        val hapGraph = dataList[0].data as HaplotypeGraph

        val dbConn = DBLoadingUtils.connection(false)
        val phgDb = PHGdbAccess(dbConn)
        readMappingDecoder = ReadMappingDecoder(phgDb)

        if (keyFile() != null) processKeyFile(hapGraph) else parentsForReadMethod(hapGraph)

        dbConn.close()

        return null
    }

    private fun processKeyFile(hapGraph: HaplotypeGraph) {
        val parentFinder = MostLikelyParents(hapGraph)
        val hapidToRefRangeMap = hapGraph.referenceRangeList().map { refrange ->
            hapGraph.nodes(refrange).map { Pair(it.id(), refrange) }
        }.flatten().toMap()

        val myWriter = PrintWriter(outputFile())
        myWriter.println("sample_name\ttotal_reads\tparent\torder\tnumber_of_reads\tcumulative_reads")
        for (sample in readKeyFile()) {
            val readMappingIdList = sample.second.map { it.toInt() }

            //create a Multimap<ReferenceRange, HapIdSetCount>
            val countsByRefrange = multimapFromMappingDataList(readMappingIdList, hapidToRefRangeMap,hapGraph)

            //call findMostLikelyParents()
            val sampleResult = parentFinder.findMostLikelyParents(countsByRefrange, maxNumberOfParents(), minCoverage())

            //write to output file
            //columns: sample_name, total_reads, parent, order, number_of_reads, cumulative_reads
            val sampleName = sample.first
            val totalReads = countsByRefrange.values().sumOf { it.count }
            val totalReadsDbl = totalReads.toDouble()
            var cumulativeSum = 0
            sampleResult.forEachIndexed { index, pair ->
                val parent = pair.first
                val order = index + 1
                val numberOfReads = pair.second
                val cumulativeReads = (pair.second + cumulativeSum) / totalReadsDbl
                cumulativeSum += pair.second
                myWriter.println("$sampleName\t$totalReads\t$parent\t$order\t$numberOfReads\t$cumulativeReads")
            }
            myWriter.flush()
        }
        myWriter.close()
    }


    private fun readKeyFile() : List<Pair<String, List<String>>> {
        val keyFile = keyFile()
        require(keyFile != null) {"Attempted to read keyFile but keyFile() was null"}
        return Files.readAllLines(Paths.get(keyFile)).drop(1).map { line ->
            val data = line.split("\t")
            Pair(data[0], data[1].split(","))
        }
    }

    private fun multimapFromMappingData(mappingId: Int,
                                        hapidToRefRangeMap: Map<Int, ReferenceRange>, graph:HaplotypeGraph):
            Multimap<ReferenceRange, HapIdSetCount> {

        val readMapping = readMappingDecoder.getDecodedReadMappingForMappingId(mappingId)

        val readMultimap = HashMultimap.create<ReferenceRange, HapIdSetCount>()
        for ( setCount in readMapping) {
            val refrange = hapidToRefRangeMap[setCount.key[0]]
            if (refrange != null) readMultimap.put(refrange, HapIdSetCount(setCount.key.toSet(), setCount.value))
        }
        return readMultimap
    }

    private fun multimapFromMappingDataList(mappingIdList: List<Int>,
                                            hapidToRefRangeMap: Map<Int, ReferenceRange>, graph:HaplotypeGraph): Multimap<ReferenceRange, HapIdSetCount> {
        require(mappingIdList.isNotEmpty()) {"Tried to decode an empty mappingDataList."}
        return if (mappingIdList.size == 1) multimapFromMappingData(mappingIdList[0], hapidToRefRangeMap,graph)
        else {
            val aggregatedReadMappingCounts = HashMultiset.create<List<Int>>()
            val refRangeToHapIdSetCounts = HashMultimap.create<ReferenceRange,HapIdSetCount>()
            for(mappingId in mappingIdList) {
                //Merge all ReadMapping counts into a single readMapping
                readMappingDecoder.getDecodedReadMappingForMappingId(mappingId).forEach { (idList, count) ->
                    aggregatedReadMappingCounts.add(idList,count)
                }
            }
            aggregatedReadMappingCounts.elementSet()
                .filter { it.isNotEmpty() }
                .forEach {
                    val referenceRange = hapidToRefRangeMap[it[0]]
                    refRangeToHapIdSetCounts.put(referenceRange,HapIdSetCount(it.toSet(),aggregatedReadMappingCounts.count(it)))
                }
            refRangeToHapIdSetCounts
        }
    }

    private fun parentsForReadMethod(hapGraph: HaplotypeGraph) {
        //same as processKeyFile except that there is always only one mapping_data per sample
        //and will need to get the sample_names and read_mapping_ids from the db
        val parentFinder = MostLikelyParents(hapGraph)
        val hapidToRefRangeMap = hapGraph.referenceRangeList().map { refrange ->
            hapGraph.nodes(refrange).map { Pair(it.id(), refrange) }
        }.flatten().toMap()

        val dbConn = DBLoadingUtils.connection(false)

        //get read_mapping_ids and line_name for the read method

        val methodQuery = "SELECT read_mapping_id, line_name, file_group_name FROM read_mapping, genotypes, methods " +
                "WHERE methods.name = '${readMethod()}' AND methods.method_id=read_mapping.method_id " +
                "AND read_mapping.genoid=genotypes.genoid"
        val mappingIdList = mutableListOf<Pair<Int,String>>()
        myLogger.info("Executing query: $methodQuery")

        dbConn.createStatement().executeQuery(methodQuery).use {
            while (it.next()) {
                mappingIdList.add(Pair(it.getInt(1), "${it.getString(2)}:${it.getString(3)}"))
            }
        }

        //TODO multithread samples (get read data from db before multithreading analysis)
        val myWriter = PrintWriter(outputFile())
        myWriter.println("sample_name\ttotal_reads\tparent\torder\tnumber_of_reads\tcumulative_reads")
        for (sample in mappingIdList) {
            val mappingIdList = listOf(sample.first)

            //create a Multimap<ReferenceRange, HapIdSetCount>
            val countsByRefrange = multimapFromMappingData(sample.first, hapidToRefRangeMap,hapGraph)

            //call findMostLikelyParents()
            val sampleResult = parentFinder.findMostLikelyParents(countsByRefrange,maxNumberOfParents(), minCoverage())

            //write to output file
            //columns: sample_name, total_reads, parent, order, number_of_reads, cumulative_reads
            val sampleName = sample.second
            val totalReads = countsByRefrange.values().sumOf { it.count }
            val totalReadsDbl = totalReads.toDouble()
            var cumulativeSum = 0
            sampleResult.forEachIndexed { index, pair ->
                val parent = pair.first
                val order = index + 1
                val numberOfReads = pair.second
                val cumulativeReads = (pair.second + cumulativeSum) / totalReadsDbl
                cumulativeSum += pair.second
                myWriter.println("$sampleName\t$totalReads\t$parent\t$order\t$numberOfReads\t$cumulativeReads")
            }
            myWriter.flush()
        }
        myWriter.close()
        dbConn.close()

    }

    override fun getIcon(): ImageIcon? {
        val imageURL = BestHaplotypePathPlugin::class.java.getResource("/net/maizegenetics/analysis/images/missing.gif")
        return if (imageURL == null) {
            null
        } else {
            ImageIcon(imageURL)
        }
    }

    override fun getButtonName(): String {
        return "LikelyParents"
    }

    override fun getToolTipText(): String {
        return "Find the most likely parents"
    }

    /**
     * The file to which the output will be written. If the
     * file exists, an error will be thrown.
     *
     * @return Out File
     */
    fun outputFile(): String {
        return outputFile.value()
    }

    /**
     * Set Out File. The file to which the output will be
     * written. If the file exists, an error will be thrown.
     *
     * @param value Out File
     *
     * @return this plugin
     */
    fun outputFile(value: String): LikelyParentsPlugin {
        outputFile = PluginParameter(outputFile, value)
        return this
    }

    /**
     * If this is set to true and the output file exists,
     * it will be overwritten without warning.
     *
     * @return Overwrite
     */
    fun overwriteOutput(): Boolean {
        return overwriteOutput.value()
    }

    /**
     * Set Overwrite. If this is set to true and the output
     * file exists, it will be overwritten without warning.
     *
     * @param value Overwrite
     *
     * @return this plugin
     */
    fun overwriteOutput(value: Boolean): LikelyParentsPlugin {
        overwriteOutput = PluginParameter(overwriteOutput, value)
        return this
    }

    /**
     * KeyFile file name.  Must be a tab separated file using
     * the following headers:
     * SampleName	ReadMappingIds
     * Additional headers will be ignored. ReadMappingIds
     * need to be comma separated for multiple values.
     * If a keyFile is not supplied then all of the read mappings
     * for the supplied method will be run.
     *
     * @return Key File
     */
    fun keyFile(): String? {
        return keyFile.value()
    }

    /**
     * Set Key File. KeyFile file name.  Must be a tab separated
     * file using the following headers:
     * SampleName	ReadMappingIds
     * Additional headers will be ignored. ReadMappingIds
     * need to be comma separated for multiple values.
     * If a keyFile is not supplied then all of the read mappings
     * for the supplied method will be run.
     *
     * @param value Key File
     *
     * @return this plugin
     */
    fun keyFile(value: String): LikelyParentsPlugin {
        keyFile = PluginParameter(keyFile, value)
        return this
    }

    /**
     * A read method. Either a read method or a keyfile containing
     * mapping ids must be supplied.
     * If both are supplied, the keyFile will be used and
     * the readMethod will be ignored.
     *
     * @return Read Method
     */
    fun readMethod(): String? {
        return readMethod.value()
    }

    /**
     * Set Read Method. A read method. Either a read method
     * or a keyfile containing mapping ids must be supplied.
     * If both are supplied, the keyFile will be used and
     * the readMethod will be ignored.
     *
     * @param value Read Method
     *
     * @return this plugin
     */
    fun readMethod(value: String): LikelyParentsPlugin {
        readMethod = PluginParameter(readMethod, value)
        return this
    }

    /**
     * The maximum number of parents to be selected.
     *
     * @return Max Parents
     */
    fun maxNumberOfParents(): Int {
        return maxNumberOfParents.value()
    }

    /**
     * Set Max Parents. The maximum number of parents to be
     * selected.
     *
     * @param value Max Parents
     *
     * @return this plugin
     */
    fun maxNumberOfParents(value: Int): LikelyParentsPlugin {
        maxNumberOfParents = PluginParameter(maxNumberOfParents, value)
        return this
    }

    /**
     * No more parents will be selected when the proportion
     * of total reads mapped to the selected parents is greater
     * than or equal to this number.
     *
     * @return Min Coverage
     */
    fun minCoverage(): Double {
        return minCoverage.value()
    }

    /**
     * Set Min Coverage. No more parents will be selected
     * when the proportion of total reads mapped to the selected
     * parents is greater than or equal to this number.
     *
     * @param value Min Coverage
     *
     * @return this plugin
     */
    fun minCoverage(value: Double): LikelyParentsPlugin {
        minCoverage = PluginParameter(minCoverage, value)
        return this
    }

}
        */