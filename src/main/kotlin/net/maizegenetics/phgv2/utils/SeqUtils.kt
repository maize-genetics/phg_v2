package net.maizegenetics.phgv2.utils

import biokotlin.seq.NucSeq
import java.io.BufferedInputStream
import java.io.ByteArrayOutputStream
import java.io.File
import org.apache.logging.log4j.LogManager

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.utils.SeqUtils")

/**
 * The file contains functions to pull data from an AGC compressed file.
 * It is assumed the AGC file is named "assemblies.agc" and is in the same folder as the tiledb datasets.
 * THis means the user only need remember a single directory, and that is sent to
 * most of the phg_v2 commands
 *
 * The user-callable commands in this file are: (Note: in the following genome = sample name)
 * 1. retrieveAgcContigs: called to pull data using the "agc getctg" command.
 *    Expected user input:  dbPath:String, ranges:List<String>
 *    The user range list should have entries that look like one of the following:
 *     contig@genome:start-end
 *     contig@genome
 *     contig:start-end (this will error if there are more than 1 genome with "contig")
 *     contig1 (this will error if there are more than 2 genome with "contig1" )
 *    Returns: Map<Pair<String,String>, NucSeq> where the key is Pair(genome, contig:start-end) or Pair(genome, contig)
 *    and the value is a Biokotlin NucSeq
 *
 * 2. retrieveAgcGenomes: called to pull data using the "agc getset" command.
 *   Expected user input:  dbPath:String, genomes:List<String>
 *   The user genome list should have sample names that are in the AGC file.
 *   Returns: Map<Pair<String,String>, NucSeq> where the key is Pair(genome, contig)
 *   and the value is a Biokotlin NucSeq
 *
 * 3. retrieveAgcData: called to pull data using    the "agc listset" or "agc listctg" command.
 *   Expected user input:  dbPath:String, agcCmdList:List<String>
 *   a.  For listset: the user agcCmdList should contain a single item: "listset"
 *       This command  returns a list of all the genomes in the AGC file.
 *   b.  For listctg: the user agcCmdList should contain "listctg" as its first entry.
 *       The remaining entries should be genome names whose contigs will be listed.
 *       If no genomes names are supplied, an error is thrown
 *   Result:  a list of strings, one per genome, in the format:
 *      genomeName:ctg1,ctg2,ctg3,ctg4
 *      e.g. LineA:1,2
 *           LineB:1,2
 *
 * 4. retrieveRefSampleName: called to get the name of the sample used as the reference in the AGC file.
 *    Expected user input: dbPath:String
 *    THe command sent to age is "agc listref <in.agc>"  THe output is written to the output stream
 *    and returned from the function as a string.
 *
 * The rest of the functions in this file are used internally by the above functions.
 *
 */

/**
 * Retrieves the name of the sample used as the reference in the AGC file.
 * @param dbPath    the folder containing assemblies.agc
 * @return the name of the sample used as the reference in the AGC file
 */
fun retrieveRefSampleName(dbPath: String,condaEnvPrefix:String):String {
    val agcFile = "${dbPath}/assemblies.agc"
    check(File(agcFile).exists()) { "File assemblies.agc does not exist in folder $dbPath- please send a valid path that indicates the parent folder for your assemblies.agc compressed file." }

    val command = if (condaEnvPrefix.isNotBlank()) arrayOf("conda","run","-p",condaEnvPrefix,"agc","listref",agcFile)
    else arrayOf("conda","run","-n","phgv2-conda","agc","listref",agcFile)

    val commandToPrint = if (command.joinToString(" ").length > 150) command.joinToString(" ").substring(0,150) else command.joinToString(" ")
    myLogger.info("retrieveRefSampleName: Running Agc Command:\n${commandToPrint}...")

    val agcProcess = ProcessBuilder(*command)
        .start()

    val agcOut = BufferedInputStream(agcProcess.inputStream, 5000000)
    val refSampleName = agcOut.bufferedReader().readLine()
    agcOut.close()
    return refSampleName
}

/**
 * Retrieves sequence from an agc data store.
 * @param dbPath    the folder containing assemblies.agc
 * @param sampleName the name of a sample in the data store
 * @param ranges    a list of pairs of Positions describing the ranges to be returned
 * @return a map of Pair<sampleName, range> to a Biokotlin NucSeq containing the sequence
 *
 *
 */
fun retrieveAgcContigs(dbPath: String, sampleName: String, ranges: List<Pair<Position,Position>>,condaEnvPrefix:String) : Map<Pair<String,String>,NucSeq> {
    val rangesString = buildRangesString(sampleName, ranges)
    return retrieveAgcContigs(dbPath,rangesString,condaEnvPrefix)

}

/**
 * Retrieves sequence from an agc data store.
 * @param dbPath    the folder containing assemblies.agc
 * @param ranges    a list of ranges to return
 * @return a map of Pair<sampleName, chrom:start-end> (or Pair<sampleName, chrom>) to a Biokotlin NucSeq containing the sequence
 *
 *    The user range list should have entries that look like one of the following:
 *     contig@genome:start-end
 *     contig@genome
 *     contig:start-end (this will error if there are more than 1 genome with "contig")
 *     contig1 (this will error if there are more than 2 genome with "contig1" )
 *
 */
fun retrieveAgcContigs(dbPath:String,ranges:List<String>,condaEnvPrefix:String):Map<Pair<String,String>,NucSeq> {

    val commands = buildAgcCommandFromList(dbPath, "getctg",ranges,condaEnvPrefix)
    val genomeChromNucSeqMap = queryAgc(commands)
    return genomeChromNucSeqMap
}

fun retrieveAgcGenomes(dbPath:String,genomes:List<String>,condaEnvPrefix:String):Map<Pair<String,String>,NucSeq> {

    val commands = buildAgcCommandFromList(dbPath, "getset",genomes,condaEnvPrefix)
    val genomeChromNucSeqMap = queryAgc(commands)
    return genomeChromNucSeqMap
}

fun buildRangesString(sampleName: String, ranges: List<Pair<Position,Position>>) : List<String> {
    val rangeStrings = mutableListOf<String>()
    for (range in ranges) {
        //ctg1@gn1:from1-to1
        //Subtracting 1 from both positions because AGC uses 0-based positions
        val rangeString = "${range.first.contig}@${sampleName}:${range.first.position-1}-${range.second.position-1}"
        rangeStrings.add(rangeString)
    }
    return rangeStrings
}

/**
 * Processes AGC commands listset or listctg
 */
fun retrieveAgcData(dbPath:String,agcCmdList:List<String>,condaEnvPrefix:String):List<String> {

    val agcFile = "${dbPath}/assemblies.agc"

    check(File(agcFile).exists()) { "File assemblies.agc does not exist in folder $dbPath- please send a valid path that indicates the parent folder for your assemblies.agc compressed file." }
    check(agcCmdList.isNotEmpty()) { "No AGC command was sent to retrieveAgcData" }

    val agcCmd = agcCmdList[0]
    val genomes = mutableListOf<String>()

    if (agcCmd == "listctg") {
        check(agcCmdList.size > 1) { "No genomes were sent to retrieveAgcData " }
        genomes.addAll(agcCmdList.subList(1,agcCmdList.size))
    }

    val command = mutableListOf<String>()
    val commandPrefix = if (condaEnvPrefix.isNotBlank()) arrayOf("conda","run","-p",condaEnvPrefix,"agc",agcCmd,agcFile)
              else arrayOf("conda","run","-n","phgv2-conda","agc",agcCmd,agcFile)
    command.addAll(commandPrefix)
    command.addAll(genomes)

    val commandToPrint = if (command.joinToString(" ").length > 150) command.joinToString(" ").substring(0,150) else command.joinToString(" ")
    myLogger.info("retrieveAgcData: Running Agc Command:\n${commandToPrint}...")
    val agcProcess = ProcessBuilder(*command.toTypedArray())
        .start()

    val agcOut = BufferedInputStream(agcProcess.inputStream, 5000000)
    try {
        // check command.  Should only be listctg or listset.
        // call appropriate function to process the output
        val queryData = when (agcCmd) {
            "listctg" -> {
                dataFromListCtgs(agcOut, genomes)
            }
            "listset" -> {
                dataFromListSet(agcOut)
            }
            else -> {
                myLogger.error("Invalid AGC command sent to retrieveAgcData: $agcCmd")
                throw IllegalStateException("Invalid AGC command sent to retrieveAgcData: $agcCmd")
            }
        }
        return queryData
    } catch (exc: Exception) {
        myLogger.error("Error reading AGC output: ${exc.message}")
        throw IllegalStateException("Error reading AGC output: ${exc.message}")
    } finally {
        agcOut.close()
    }
}

/**
 * The dataFromListSet method parses the AGC output from a listset command.
 */
fun dataFromListSet(agcOut:BufferedInputStream):List<String>{
    val data = mutableListOf<String>()  // this is the list that will be returned

    try {
        agcOut.bufferedReader().use { br ->
            var line = br.readLine()
            while (line != null) {
                // AGC often returns blank line at the end - skip it.
                if (line != "") data.add(line)
                line = br.readLine()
            }
        }
    } catch (exc: Exception) {
        myLogger.error("Error reading AGC output: ${exc.message}")
    } finally {
        agcOut.close()
    }
    return data
}

/**
 * The dataFromListCtgs method parses the AGC output from a listctg command.
 */
 fun dataFromListCtgs(agcOut:BufferedInputStream, genomes:List<String>):List<String> {
     val dataMap = mutableMapOf<String, MutableList<String>>()
     try {
         agcOut.bufferedReader().use { br ->
             var line = br.readLine()
             var currentGenome = ""
             while (line != null) {
                 // if line matches one of the genomes in the genomes list, create a new list for that genome
                 // and add the line to that list
                 if (line == "") {
                     // AGC has blank lines at times - skip them
                     line = br.readLine()
                     continue
                 }
                 line = line.trim()
                 if (genomes.contains(line)) {
                     currentGenome = line
                 } else {
                     // if the line does not match one of the genomes, add it to the list for the last genome
                     var dataList = dataMap[currentGenome]
                     if (dataList == null) {
                         dataList = mutableListOf<String>()
                         dataMap[currentGenome] = dataList
                     }
                     dataList.add(line)
                 }
                 line = br.readLine()
             }
         }
     } catch (exc: Exception) {
         myLogger.error("Error reading AGC output: ${exc.message}")
     } finally {
         agcOut.close()
     }
     // Turn the datamap into a List<String> where each string is composed of the map key followed by a colon and then
     // a comma separated list of the contigs for that genome
     val dataList = mutableListOf<String>()
     for (genome in dataMap.keys) {
         val contigList = dataMap[genome]
         val contigString = contigList?.joinToString(",")
         dataList.add("$genome:$contigString")
     }
     return dataList
 }


/**
 * This builds the command that is sent to AGC.  It creates the initial
 *  command, which sets the conda environment, then adds the AGC file
 *  path, then follows with the commands the user has put together in a list
 *
 *  The agcCmd parameter is the specific agc command to run, e.g. "getctg",  or "getset"
 *  This method is called from retrieveAgcContigs() and retrieveAgcGenomes()
 */
fun buildAgcCommandFromList(dbPath:String, agcCmd:String, ranges:List<String>, condaEnvPrefix:String): Array<String> {

    // Validate the dbPath. We check if the folder exists, and if it contains the file assemblies.agc
    val agcFile = "${dbPath}/assemblies.agc"

    check(File(agcFile).exists()) { "File assemblies.agc does not exist in folder $dbPath- please send a valid path that indicates the parent folder for your assemblies.agc compressed file." }

    val command = mutableListOf<String>()
    val commandPrefix = if (condaEnvPrefix.isNotBlank()) arrayOf("conda","run","-p","phgv2-conda") else arrayOf("conda","run","-n","phgv2-conda")
    val commandData = arrayOf("agc",agcCmd,agcFile)
    command.addAll(commandPrefix)
    command.addAll(commandData)
    command.addAll(ranges)

    return command.toTypedArray()
}

/**
 * Retrieve AGC data from a list of commands sent via ProcessBuilder
 * This function should only be called when the AGC command is "getctg" or "getset"
 * as it specifically returns a Map<String, NucSeq>.  Other AGC commands processed
 * in this repository will return a List<String> and are processed in retrieveAgcData
 *
 * What AGC does is take the entire idline and append the range to it, (excluding the sample name).
 * When the fasta looks like:
 *    >1
 * And the query is:
 *     1@LineA:1-1000
 * AGC will return:
 *     >1:1-1000
 * When we have the sample name including in the id line, e.g
 *    >1 sampleName=LineA
 *  And the query is:
 *    1@LineA:1-1000
 *  AGC will return:
 *    >1 sampleName=LineA:1-1000
 *  So we have to look for ":" or end of line to get th end of the sample name
 */
fun queryAgc(commands:Array<String>):Map<Pair<String,String>,NucSeq> {
    check(commands.size > 0) { "Error:  No commands sent to queryAgc!" }
    val commandToPrint = if (commands.joinToString(" ").length > 150) commands.joinToString(" ").substring(0,150) else commands.joinToString(" ")
    myLogger.info("queryAgc: Running Agc Command:\n${commandToPrint}...")

    // .redirectError(ProcessBuilder.Redirect.INHERIT) redirects to stdout
    // For now, running without redirecting output or error.  The errorStream will be
    // checked and handled at the end.
    val agcProcess = ProcessBuilder(*commands)
        .start()

    val agcOut = BufferedInputStream(agcProcess.inputStream, 5000000)

    // Note AGC does not include the genome in the idline, so we
    // can end up with the same idline multiple times, e.g. if we want chrom 1 from
    // multiple genomes, all the idlines will be ">1" and we will only get the last one
    // when processing into a map.
    // It is the same if we give a region for the chrom, e.g. 1@LineA:1-1000 and 1@LineB:1-1000,
    // we will get the same idline multiple times.
    val genomeChromNucSeq = HashMap<Pair<String,String>,NucSeq>()
    try {
        agcOut.bufferedReader().use { br ->
            var currChrom = "-1"
            var currSeq = ByteArrayOutputStream()
            var currGenome = "-1"
            var currRange = ""
            var line = br.readLine()
            while (line != null) {
                line = line.trim()
                if (line.startsWith(">")) {
                    if (currChrom != "-1") {
                        // finished with this chromosome's sequence
                        //myLogger.info("queryAgc: finished chrom $currChrom for genome $currGenome")
                        val chromPlusRange = if (currRange == "") currChrom else "${currChrom}:${currRange}"
                        genomeChromNucSeq.put(Pair(currGenome,chromPlusRange),NucSeq(currSeq.toString()))
                    }
                    // reset chromosome name and sequence, begin processing next chrom
                    // Split on ALL whitespace, to handle tabs etc between the chrom name and other data on the line
                    currChrom = line.substring(1).split("\\s+".toRegex())[0]

                    val genomeStart = line.indexOf("sampleName=")
                    if (genomeStart < 0) {
                        // This is an error as we need the sampleName to map the genome to the
                        // samplename/genome in the VCF files.
                        myLogger.error("queryAgc: Error:  No sampleName= found in idline: returned from AGC query.")
                        myLogger.error("Please run the prepare-assemblies command to add sampleName to your fasta files.")
                        myLogger.error("Then delete the assemblies.agc file from you tiledb folder and recreate the agc compressed file via the agc-compress command using the new prepared fastas.")
                        myLogger.error("This is necessary to associate the VCF samples with the AGC compressed genomes.")
                        throw IllegalStateException("Error:  No sampleName found in idline: $line")
                    }
                    // indexof gives the index from the start of the string, but it starts looking
                    // from where you tell it to look, in this case, from "genomeStart"
                    val sampleNameEnd = if (line.indexOf(":", genomeStart) > 0) line.indexOf(":", genomeStart) else line.length

                    currRange = if (sampleNameEnd == line.length ) "" else line.substring(sampleNameEnd+1)
                    // genomeStart is the start of "sampleName="  We add 11 to get past this to the actual samplename
                    currGenome = line.substring(genomeStart+11, sampleNameEnd)

                    currSeq = ByteArrayOutputStream()
                } else {
                    currSeq.write(line.toByteArray())
                }
                line = br.readLine()
            }
            if (currSeq.size() > 0) {
                myLogger.info("queryAgc: finished chrom $currChrom")
                val chromPlusRange = if (currRange == "") currChrom else "${currChrom}:${currRange}"
                genomeChromNucSeq.put(Pair(currGenome,chromPlusRange),NucSeq(currSeq.toString()))
            }
        }
    } catch (exc:Exception) {
        myLogger.error("Error:  Exception in queryAgc: ${exc.message}")
        throw IllegalStateException("Error:  Exception in queryAgc: ${exc.message}")
    }

    // Check for any errors
    val agcError = BufferedInputStream(agcProcess.errorStream, 5000000)
    val errors = agcError.bufferedReader().readLines()
    if (errors != null && errors.size > 0) {
        myLogger.error("queryAgc: errors found in errorStream: ${errors.size}")
        throw IllegalArgumentException("Error running AGC command: ${commands.joinToString(" ")}\nError: $errors")
    }
    return genomeChromNucSeq
}
