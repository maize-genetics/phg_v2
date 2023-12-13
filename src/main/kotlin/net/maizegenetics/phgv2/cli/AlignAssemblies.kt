package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import org.apache.logging.log4j.LogManager
import java.io.File
import java.lang.management.ManagementFactory
import javax.management.MBeanServer
import javax.management.ObjectName

/**
 * This will align assemblies to a reference genome.
 * Uses anchorwave's proali to align, which handles
 * genome alignment with relocation variation,
 * chromosome fusion or whole genome duplication.
 *
 * This function allows for multiple assembly alignments to be run
 * in parallel.  Users may specify number of alignments to run in parallel
 * and the total number of threads available to be split between the alignments.
 * Anchorwave takes 10G of RAM for the dynamic program aspects, and then another
 * 10G+ per thread, depending on the processor.
 *
 * The table below shows the memory usage for a single assembly alignment
 * based on processor type: (From Baoxing Song, 2023-11-10)
 *
 *  ---------------------------------------------
 *  | Processor | peak memory (Gb) | wall time  |
 *  ---------------------------------------------
 *  |  SSE2    |    20.1          | 26:47:17    |
 *  ---------------------------------------------
 *  |  SSE4.1  |    20.6          | 24:05:07    |
 *  ---------------------------------------------
 *  |  AVX2    |    20.1          | 21:40:00    |
 *  ---------------------------------------------
 *  |  AVX512  |    20.1          | 18:31:39    |
 *  ---------------------------------------------
 *  | ARM      |    34.2          | 18:08:57    |
 *  ---------------------------------------------
 *
 *  For a machine with 128G RAM, you would use 10 threads
 *  divided between your parallel runs.  For example, if you
 *  have 10 assemblies to align, you could run 2 parallel
 *  alignments, each using 5 threads.
 *
 */
class AlignAssemblies : CliktCommand(help="Align assemblies using anchorwave") {

    private val myLogger = LogManager.getLogger(AlignAssemblies::class.java)

    val gff by option(help = "Full path to the reference gff file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--gff must not be blank"
            }
        }

    val referenceFile by option(help = "Full path to reference fasta file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }

    val assemblies by option(
        "-a",
        "--assemblies",
        help = "File containing list of assemblies to align, 1 per line, full path to file"
    )
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--assemblies must not be blank"
            }
        }

    val outputDir by option("-o", "--output-dir", help = "Directory where temporary and final files will be written")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir must not be blank"
            }
        }

    val totalThreads by option(help = "Number of threads available.  These will be split among the alignments that are run in parallel")
        .int()
        .default(0)

    val inParallel by option(
        help = "Number of assemblies to simultaneously process. " +
                "If you have 10 threads and the in-parallel value is 2, then 2 assemblies will be aligned at a time, each using 5 threads." +
                "The anchorwave application can take up to 30G per thread for each assembly processed, plus some overhead." +
                "Consider this memory factor when providing values for the total-threads and in-parallel."
    )
        .int()
        .default(0)

    val refMaxAlignCov by option(help = "Anchorwave proali parameter R, indicating reference genome maximum alignment coverage.")
        .int()
        .default(1)

    val queryMaxAlignCov by option(help = "Anchorwave proali parameter Q, indicating query genome maximum alignment coverage.")
        .int()
        .default(1)

    data class InputChannelData(
        val refFasta: String,
        val asmFasta: String,
        val outputDir: String,
        val gffFile: String,
        val refSamOutFile: String,
        val maxRuns:Int,
        val threadsPerRun:Int
    )

    override fun run() {

        // Returns Pair<Int, Int> where the first value is the number of parallel alignments, the second is threadsPerAlignment
        val runsAndThreads = calculatedNumThreadsAndRuns(totalThreads, inParallel, assemblies)

        // create CDS fasta from reference and gff3 file
        val cdsFasta = "$outputDir/ref.cds.fasta"

        createCDSfromRefData(referenceFile, gff, cdsFasta, outputDir)

        // create list of assemblies to align from the assemblies file
        val assembliesList = File(assemblies).readLines().filter { it.isNotBlank() }

        // run minimap2 for ref to refcds
        val justNameRef = File(referenceFile).nameWithoutExtension
        val samOutFile = "${justNameRef}.sam"
        val refSamOutFile = "${outputDir}/${samOutFile}"

        // For minimap2, we will use the number of processors available as the number of threads.
        val builder = ProcessBuilder(
            "conda", "run", "-n", "phgv2-conda", "minimap2", "-x", "splice", "-t", runsAndThreads.second.toString(), "-k", "12",
            "-a", "-p", "0.4", "-N20", referenceFile, cdsFasta, "-o", refSamOutFile
        )
        val redirectError = "$outputDir/minimap2Ref_error.log"
        val redirectOutput = "$outputDir/minimap2Ref_output.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("Ref minimap Command: " + builder.command().joinToString(" "));

        val process = builder.start()
        val error = process.waitFor()
        if (error != 0) {
            myLogger.error("minimap2 for $referenceFile run via ProcessBuilder returned error code $error")
            throw IllegalStateException("Error running minimap2 for reference: $error")
        }

        runAnchorWaveMultiThread(referenceFile, assembliesList, cdsFasta, gff, refSamOutFile,runsAndThreads)

    }

    /**
     * This method uses a java managment facotry bean to determine the
     * amount of system memory on the hosting machine.
     */
    fun getSystemMemory():Double {
        val mBeanServer: MBeanServer = ManagementFactory.getPlatformMBeanServer()
        // attribute is in bytes
        val attribute =
            mBeanServer.getAttribute(ObjectName("java.lang", "type", "OperatingSystem"), "TotalPhysicalMemorySize")
        // translate bytes to gigabutes
        val memoryGibi = attribute as Long / 1e9
        // and now to gibibytes, as anchorwave uses GiB
        val memoryGigi = memoryGibi * 0.93
        myLogger.info("getSystemMemory: Total system memory: $attribute Bytes, $memoryGibi GB, $memoryGigi GiB")
        return memoryGigi
    }

    /**
     * This function calculates the number of available processors, and the maxium
     * system memory avilable for this machine.  It returns a Pair<Int, Int> where
     * the first value is the number of alignments to run in parallel, and the second
     * value is the number of threads to use for each alignment.
     *
     * An algorithm is used that selects a middle value between the number of threads
     * and the number of alignments to run in parallel.
     *
     */
    fun calculatedNumThreadsAndRuns(totalThreads:Int, inParallel:Int, assemblies:String): Pair<Int, Int> {
        // If totalThreads or inParallel are 0, it means the user did not specify them
        // In that case we calculate these values based on the number of processors available
        // and the amount of free memory available.
        val processors = Runtime.getRuntime().availableProcessors() - 2 // leave 2 processors for the OS
        val systemMemory = getSystemMemory()

        myLogger.info("calculateNumThreadsAndRuns: systemMemory: $systemMemory, processors: $processors")

        // If the user did not specify the number of threads to use, we will use all available
        // processors.  If the user did specify the number of threads, we will use that number
        // of threads, but we will not exceed the number of processors available.
        val totalThreadsToUse = if (totalThreads > 0) {
            if (totalThreads > processors) {
                myLogger.warn("The number of threads specified ($totalThreads) exceeds the number of processors available ($processors).  Using $processors threads.")
                processors
            } else {
                totalThreads
            }
        } else {
            processors
        }

        myLogger.info("calculateNumThreadsAndRuns: totalThreadsToUse: $totalThreadsToUse")
        // Per Baoxing's chart, it takes just over 20G/thread to run anchorwave on a sample maize genome.
        // THe number of threads that can be run simultaneously is limited by the amount of
        // free memory available.  We will use 21G/thread as the memory requirement, and
        // calculate the number of threads that can be run simultaneously based on the
        // amount of system memory available. This calculation assumes the user has access to
        // all memory and processors on the machine.

        val concurrentThreads = (systemMemory/21).toInt()
        myLogger.info("calculateNumThreadsAndRuns: max concurrent threads: $concurrentThreads")
        if (concurrentThreads < 1) {
            //TODO: Lcj - do we still need this code?
            // systemMemory may be low when running CI or junit tests.
            // The test files are also very small.
            // This clause allows junits to run.
            myLogger.warn("There is not enough systemMemory to run anchorwave.  Free memory: $systemMemory")
            myLogger.warn("will attempt to run one alignment at a time, using just a single thread.")
            return(Pair(1,1))
        }

        // Calculate how many threads we can run at the same time (based on memory availability)
        val totalConcurrentThreads = if (concurrentThreads > totalThreadsToUse) {
            totalThreadsToUse
        } else {
            concurrentThreads
        }

        // Now that we know how many threads can be run concurrently, we need to
        // determine how many parallel alignments to run, and how many threads each
        // one gets.  If the user specified the number of parallel alignments to run,
        // we will use that number, or the number of assemblies in the list, whichever
        // is smaller.

        val numAssemblies = File(assemblies).readLines().filter { it.isNotBlank() }.size
        // This needs to return a Pair<Int, Int> where the first value is the number of alignments, the seconds is threadsPerAlignment
        val runsAndThreads = if (inParallel > 0) {
            if (inParallel > totalConcurrentThreads) {

                myLogger.warn("The number of parallel alignments specified ($inParallel) exceeds the number of threads that can be run concurrently ($totalConcurrentThreads).  Will run $totalConcurrentThreads concurrent alignments with 1 thread each.")
                if (numAssemblies < totalConcurrentThreads){
                    // there are fewer assemblies than threads available, so split the number of threads
                    // evenly among the assemblies.  This may leave some threads unused.
                    Pair(numAssemblies, totalConcurrentThreads/numAssemblies)
                } else {
                    // we have as many concurrent aligments as there are threads available, so each alignment gets 1 thread
                    Pair(totalConcurrentThreads, 1)
                }

            } else {
                // InParallel was defined by the user and is <= totalConcurrentThreads
                // THis is adjusted lower if there are fewer assemblies to align than the number of parallel alignments
                // Number of concurrent alignments is "inParallel", and number of threads-per-alignment becomes
                //  totalConcurrentThreads/inParallel
                val concurrentAlignments = if (numAssemblies < inParallel) numAssemblies else inParallel
                Pair(concurrentAlignments, totalConcurrentThreads / concurrentAlignments)
            }
        } else {
            // User did not specify the number of alignments they want run in parallel (the inParallel param)
            // Need to do calculations here to determine the best values for
            // concurrent alignments and number of threads per alignment
            maximizeRunsAndThreads(totalConcurrentThreads, numAssemblies)
        }
        myLogger.info("calculatedNumThreadsAndRuns: returning runsAndThreads values: $runsAndThreads")
        return runsAndThreads
    }

    /**
     * This considers the number of threads on the machine, the number of threads user may have
     * specified, and the number of assemblies to align.  It calculates the a middle value balancing
     * number of threads and number of concurrent alignments.
     *
     * If the
     */
    fun maximizeRunsAndThreads(totalConcurrentThreads:Int, totalAssemblies:Int): Pair<Int, Int> {
        // This returns a Pair<Int, Int> where the first value is the number of runs, the seconds is threadsPerRun
        // THe maximum number of current runs is the total number of threads available,
        // which would be each assembly getting a single thread.
        // We believe it is more efficient to run fewer assemblies at a time, each with more threads.
        // What this code does is calculate the middle of this.  For example: If there are 20 assemblies
        // and 10 threads, our options are:
        // 1 run: 10 threads
        // 2 runs: 5 threads
        // 3 runs: 3 threads
        // 4 runs: 2 threads
        // 5 runs: 2 threads
        // 6 or higher runs: 1 thread
        // This code will return 3 runs, 3 threads each.
        myLogger.info("maximizeRunsAndThreads: totalConcurrentThreads: $totalConcurrentThreads, totalAssemblies: $totalAssemblies")
        val assembliesToThreads = mutableMapOf<Int, Int>()

        // This loop says if each assembly gets "numThreads", how many concurrent runs can we do?
        for (numThreads in 1..totalConcurrentThreads) {

            val numRuns =totalConcurrentThreads/numThreads
            if (numRuns > totalAssemblies) continue // invalid, number of runs is greater than the number of assemblies
            val currentThreads = assembliesToThreads[numRuns]
            // if currentThreads is not null and is > than numThreads, ignore.
            // otherwise, replace this entry
            if (currentThreads == null || currentThreads < numThreads) {
                assembliesToThreads[numRuns] = numThreads
            }
        }
        // we should now have a map with the highest number of threads for each number of runs
        //At this point, we pick
        // 1.  if only 1 entry, use that
        // 2.  if are only 2 entries, use the one with the highest number of threads
        // 3.  if there are > 3 entries, drop the one with the lowest number of runs and the one with the highest number of runs.
        // Repeat until there are 2 entries or fewer entries left.


        myLogger.info("maximizeRunsAndThreads: potential run/thread combinations:")
        myLogger.info("numAlignments\tthreadsPerAlignments")
        assembliesToThreads.forEach { (numRuns, threadsPerRun) ->
            myLogger.info("$numRuns\t$threadsPerRun")
        }

        // 1.  if only 1 entry, use that
        if (assembliesToThreads.size == 1) {
            val entry = assembliesToThreads.entries.first()
            myLogger.info("Running ${entry!!.key} runs with ${entry.value} threads per runs")
            return Pair(entry.key, entry.value)
        } else if (assembliesToThreads.size == 2) {
            // 2.  if are only 2 entries, use the one with the highest number of concurrent assemblies
            val entry = assembliesToThreads.entries.maxByOrNull { it.key }
            myLogger.info("Running ${entry!!.key} runs with ${entry.value} threads per runs")
            return Pair(entry.key, entry.value)
        } else {
            // 3.  if there are > 3 entries, drop the one with the lowest number of runs and the one with the highest number of runs.
            // Repeat then there are 2 entries or fewer entries left.
            while (assembliesToThreads.size > 2) {
                val minEntry = assembliesToThreads.entries.minByOrNull { it.key }
                val maxEntry = assembliesToThreads.entries.maxByOrNull { it.key }
                assembliesToThreads.remove(minEntry!!.key)
                assembliesToThreads.remove(maxEntry!!.key)
            }
            // 2.  if are only 2 entries, use the one with the highest number of concurrent alignments
            // This may mean we are not using all the threads available, but it avoids the case of choosing
            // to run 1 alignment with all threads, rather than 2 alignments with hald the threads when
            // these are our only choices.

            val entry = assembliesToThreads.entries.maxByOrNull { it.key }
            myLogger.info("Running ${entry!!.key} concurrent alignments with  ${entry!!.value} threads for each run")
            return Pair(entry.key, entry.value)
        }
    }

    private fun createCDSfromRefData(refFasta: String, gffFile: String, cdsFasta: String, outputDir: String): Boolean {

        // val command = "anchorwave gff2seq -r ${refFasta} -i ${gffFile} -o ${cdsFasta} "
        // Need to set the conda environment here to access anchorwave
        val builder = ProcessBuilder(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "anchorwave",
            "gff2seq",
            "-r",
            refFasta,
            "-i",
            gffFile,
            "-o",
            cdsFasta
        )

        val redirectOutput = "$outputDir/anchorwave_gff2seq_output.log"
        val redirectError = "$outputDir/anchorwave_gff2seq_error.log"
        builder.redirectOutput(File(redirectOutput))
        builder.redirectError(File(redirectError))

        myLogger.info("createCDSfromRefData command:" + builder.command().joinToString(" "))
        try {
            val process = builder.start()
            val error = process.waitFor()
            if (error != 0) {
                myLogger.error("createCDSfromRefData run via ProcessBuilder returned error code $error")
                return false
            }
            return true
        } catch (e: Exception) {
            myLogger.error("Error: could not execute anchorwave command. Run anchorwave manually and retry.")
            return false
        }

    }

    /**
     * This will run anchorwave for each assembly in the assemblies list.
     * Each anchorwave process will be run in a separate thread.
     */
    private fun runAnchorWaveMultiThread(
        refFasta: String,
        assemblies: List<String>,
        cdsFasta: String,
        gffFile: String,
        refSamOutFile: String,
        runsAndThreads: Pair<Int, Int> // Pair(concurrentRuns, threadsPerRun)
    ) {
        runBlocking {
            // Setup
            val inputChannel = Channel<InputChannelData>()

            // The input channel gets data needed to run minimap2 and align with anchorwave
            launch {
                myLogger.info("Adding entries to the inputChannel:")
                assemblies.forEach { asmFile ->

                    // Column names were checked for validity above
                    myLogger.info("Adding: $asmFile for processing")
                    inputChannel.send(InputChannelData(refFasta, asmFile, outputDir, gffFile, refSamOutFile, runsAndThreads.first, runsAndThreads.second))
                }
                myLogger.info("Done Adding data to the inputChannel:")
                inputChannel.close() // Need to close this here to show the workers that it is done adding more data
            }

            // Do not need a coroutine that "joins" the threads as they will all
            // terminate above when there is no more data on the input channel

            // This calls anchorwave's proali, and minimap2 scripts to process the alignments
            // The number of worker threads should equal the number of concurrent runs
            // Inside alignAssembly(), each alignment will get threadsAndRuns.second, which is the number
            // of threads allocated for each alignment.  This data is included in the InputChannelData record
            val workerThreads = (1..runsAndThreads.first).map { run ->
                launch { alignAssembly(inputChannel, cdsFasta, gffFile) }
            }

        }
    }

    private suspend fun alignAssembly(inputChannel: Channel<InputChannelData>, cdsFasta: String, gffFile: String) =
        withContext(
            Dispatchers.Default
        ) {
            //val threadsPerRun = totalThreads / inParallel
            //println("alignAssembly: totalThreads: $totalThreads, inParallel: $inParallel, threadsPerRun: $threadsPerRun")

            for (assemblyEntry in inputChannel) {
                // Column names were checked for validity above
                val justName = File(assemblyEntry.asmFasta).nameWithoutExtension
                val samShort = "${justName}.sam"
                val asmSamFile = "${assemblyEntry.outputDir}/${samShort}"

                myLogger.info("alignAssembly: asmFileFull: ${assemblyEntry.asmFasta}, outputFile: $asmSamFile , threadsPerRun: ${assemblyEntry.threadsPerRun}")

                val builder = ProcessBuilder(
                    "conda",
                    "run",
                    "-n",
                    "phgv2-conda",
                    "minimap2",
                    "-x",
                    "splice",
                    "-t",
                    assemblyEntry.threadsPerRun.toString(),
                    "-k",
                    "12",
                    "-a",
                    "-p",
                    "0.4",
                    "-N20",
                    assemblyEntry.asmFasta,
                    cdsFasta,
                    "-o",
                    asmSamFile
                )

                val redirectError = "${assemblyEntry.outputDir}/minimap2_${justName}_error.log"
                val redirectOutput = "${assemblyEntry.outputDir}/minimap2_${justName}_output.log"
                myLogger.info("redirectError: $redirectError")
                builder.redirectOutput(File(redirectOutput))
                builder.redirectError(File(redirectError))
                myLogger.info(" begin minimap assembly Command: " + builder.command().joinToString(" "))

                val process = builder.start()
                val error = process.waitFor()
                if (error != 0) {
                    myLogger.error("minimap2 for assembly ${assemblyEntry.asmFasta} run via ProcessBuilder returned error code $error")
                    throw IllegalStateException("alignAssembly: error running minimap2 for ${justName}: $error")
                }
                // We have the SAM File, call proali to align with anchorwave
                runAnchorwaveProali(
                    gffFile,
                    assemblyEntry.refFasta,
                    assemblyEntry.asmFasta,
                    cdsFasta,
                    assemblyEntry.refSamOutFile,
                    asmSamFile,
                    assemblyEntry.threadsPerRun.toString()
                )
            }
        }

    private fun runAnchorwaveProali(
        gffFile: String,
        refFasta: String,
        asmFasta: String,
        cdsFasta: String,
        refSam: String,
        asmSam: String,
        threadsPerRun: String
    ) {

        val justNameAsm = File(asmFasta).nameWithoutExtension
        val justNameRef = File(refFasta).nameWithoutExtension

        val anchorsproFile = "${outputDir}/${justNameAsm}_${justNameRef}.anchorspro"
        // Using just the assembly name.  This facilitates programmatically creating
        // a GVCF keyfile from the maf keyfiles.  It will be understood that the maf
        // file name is <assemblyFastaNoExtension>.maf
        val outputFile = "${outputDir}/${justNameAsm}.maf"
        val builder = ProcessBuilder(
            "conda",
            "run",
            "-n",
            "phgv2-conda",
            "anchorwave",
            "proali",
            "-i",
            gffFile,
            "-r",
            refFasta,
            "-as",
            cdsFasta,
            "-a",
            asmSam,
            "-ar",
            refSam,
            "-s",
            asmFasta,
            "-n",
            anchorsproFile,
            "-R",
            refMaxAlignCov.toString(),
            "-Q",
            queryMaxAlignCov.toString(),
            "-t",
            threadsPerRun,
            "-o",
            outputFile
        )

        val redirectError = "${outputDir}/proali_${justNameAsm}_outputAndError.log"
        myLogger.info("redirectError: $redirectError")

        // NOTE: anchowave proali has an output file parameter, unlike minimap2, which uses
        // command line redirection (ie >) to write the .SAM file
        //
        // Do not use the anchorwave output file parameter as the message redirect output
        // for ProcessBuilder  or you'll end up with the text "AnchorWave done!" in the middle of
        // a reference sequence, causing problems in parsing, and causing the sequence itself to be
        // on a separate line from the data describing it.  Instead, send all program output messages
        // to the error file.
        // Here, we  write all program message output (ie non MAF file output) to a single file
        // named as above.

        builder.redirectOutput(File(redirectError))
        builder.redirectError(File(redirectError))
        myLogger.info(
            "runAnchorwaveProali proali Command for ${justNameAsm}: " + builder.command().joinToString (" ")
        )

        try {
            val process = builder.start()
            val error = process.waitFor()
            if (error != 0) {
                myLogger.error("proali for assembly $asmFasta run via ProcessBuilder returned error code $error")
                throw IllegalStateException("runAnchorwaveProali: error running proali for $justNameAsm")
            }
        } catch (e: Exception) {
            myLogger.error("Error: could not execute anchorwave command. Run anchorwave manually and retry.")
        }
    }

}