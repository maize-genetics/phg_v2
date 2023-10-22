package net.maizegenetics.phgv2.cli

import biokotlin.seq.NucSeq
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.groups.OptionGroup
import com.github.ajalt.clikt.parameters.groups.cooccurring
import com.github.ajalt.clikt.parameters.groups.provideDelegate
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.options.validate
import com.github.ajalt.clikt.parameters.types.int
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream
import java.io.BufferedReader
import java.io.ByteArrayOutputStream

/**
 * generic command that pulls data from an AGC compressed file
 * User may specify a list of sample names, list of contigs and
 * a range of positions to pull from the file.
 *
 * SampleName is the only required parameter.
 *
 * If samplename is not specified, all will be pulled
 * If contig is not specified, all will be pulled
 *
 * If range IS specified, contig and sample must be specified.
 *
 * The output for now will be streamed back to user.
 */


class AgcGet : CliktCommand(help="Pull data from an AGC compressed file") {

    private val myLogger = LogManager.getLogger(AgcGet::class.java)

    val dbPath by option(help = "Folder name where AGC compressed filelives.  Should be the same parent folder where tiledb datasets are stored.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--db-path must not be blank"
            }
        }
    val sampleNames by option (help = "Comma separated list of samples from which data will be pulled.  If not specified, all samples will be pulled.")
        .default("")
        .validate{
            require(it.isNotBlank()) {
                "--sample-names must not be blank"
            }
        }

    val contigs by option (help = "Comma separated list of contig names to pull.  If not specified, all contigs will be pulled.")
        .default("")

    val start by option(help = "Start 1-based position of range from which data will be pulled.  If specified, end must also be specified." )
        .int()
        .default(0).validate {
            require(it >= 0 ) {
                "if --end is present, start must be greater than 0 and less than or equal to end"
            }
        }

    val end by option (help = "End 1-based position of range from which data will be pulled.  If specified, start must also be specified.")
        .int()
        .default(0).validate{
            require(it >= 0 && it >= start) {
                "if --start is present, end must be greater than or equal to --start"
            }
        }


    override fun run() {
        println("AgcGet - in run !!")
        validateParameters()
        val command = buildAgcCommandFromParams(dbPath, sampleNames, contigs, start, end)
        println("command = ${command}")
    }

    private fun validateParameters() {
        // Start and end were given defaults of 0 if not passed by user
        check(start <= end) {"Error, start value ${start} is greater than end value ${end}"}

        // Verify if range is present, then contig is present
        if (start > 0 && contigs=="") {
            throw IllegalArgumentException("Contig must be specified when start/end positions are given")
        }

        if (end > 0 && start < 1) {
            // This was hard to catch with parameter checking.
            // we need a start value if we have an end value, and as these are 1-based,
            // the value must be > 0.  End > 0 is checked on the parameter.
            throw IllegalArgumentException("Start parameter must be specified when end parameter is specified and it must be 1 or greater")
        }
    }

    // This function builds the agc command from the user supplied
    // parameters.  It is meant to be called from this class's run() program
    fun buildAgcCommandFromParams(dbPath:String, sampleNames:String, contigs:String, start:Int, end:Int ):Array<String> {

        // Set the beginning data for the ProcessBuilder command
        val command = mutableListOf<String>()
        val commandPrefix = arrayOf("conda","run","-n","phgv2-conda","agc","getctg",dbPath)
        command.addAll(commandPrefix)
        var sampleList = sampleNames.split(",").map { it.trim() }
        if (contigs == "") {
            // no contigs, so no positions.  user just wants the genomes in a file
            // split the sampleNames(genomes) into a list and add to the command.
            // this creates a  command of the format: (out.fa is the redirect for ProcessBuilder
            //     agc getctg dbpath sampleName1 sampleName2 ... > out.fa
            command.addAll(sampleList)
        } else if (start > 0) {
            // we have contigs and positions. If we have a range, we apply
            // the same range to all contigs.  If user wanted different ranges
            // for different contigs, they would need to callbuildAgcCommandFromList()
            // from some other program
            // It is not expected there would be multiple contigs and multiple
            // sampleNames with the same ranges.  But if there are, all contigs
            // and ranges are applied to all sampleNames
            val contigList = contigs. split(",").map { it.trim() }
            sampleList.forEach { sample ->
                contigList.forEach{ contig ->
                    val data = "${contig}@${sample}:${start}-${end}"
                    command.add(data)
                }
            }
        } else {
            // This should be the case with genomes (sampleNames) and contigs,
            // but no ranges specified.
            val contigList = contigs. split(",").map { it.trim() }
            sampleList.forEach { sample ->
                contigList.forEach{ contig ->
                    val data = "${contig}@${sample}"
                    command.add(data)
                }
            }
        }

        // Command list has been created, it is ready to be loaded to agc.
        return command.toTypedArray()
        //return retrieveAgcData(command.toTypedArray())
    }

    // This is available for calling from other programs.  It will return an array
    // of commands that can be sent to the agc here.  It's intended use is for user
    // who need to have different ranges pulled from different genome/contigs,
    // The user input would be a list of something like:
    //    ctg1@gn1:from1-to1
    //    ctg2@gn2:from2-to2
    //    etc.
    fun buildAgcCommandFromList(dbPath:String, ranges:List<String>): Array<String> {
        // this is meant to be called from other programs.  it will return an array
        // of commands that can be sent to agc
        val commands = mutableListOf<String>()

        // TODO: finish this !!
        // Perhaps this  calls agc to load
        // and returns that data
        return commands.toTypedArray()
    }

    // function to retrieve the data from agc
    //fun retrieveAgcData(commands:Array<String>):BufferedInputStream {
    fun retrieveAgcData(commands:Array<String>):Map<String, NucSeq> {
        check(commands.size > 0) { "Error:  No commands sent to retrieveAgcData!" }
        myLogger.info("Running Agc Command:\n${commands.joinToString(" ")}")
        val agcProcess = ProcessBuilder(*commands)
            .redirectError(ProcessBuilder.Redirect.INHERIT)
            .start()
        // Do we return the bufferedInputStream below, or process into a map of chrom to NucSeq?
        val agcOut = BufferedInputStream(agcProcess.inputStream, 5000000)

        // Or ... do we want to process the results into a Map<Idline, NucSeq>
        // THe problem is, AGC does not include the genome in the idline, so we
        // can end up with the same idline multiple times, e.g. if we want chrom 1 from
        // multiple genomes, all the idlines will be ">1" and we will only get the last one
        // when processing into a map.
        // It is the same if we give a region for the chrom, e.g. 1:1-1000, we will get
        // the same idline multiple times.
        val chromNucSeqMap = HashMap<String, NucSeq>()
        try {
            agcOut.bufferedReader().use { br ->
                var currChrom: String = "-1"
                var currSeq = ByteArrayOutputStream()
                var line = br.readLine()
                while (line != null) {

                    line = line.trim()
                    if (line.startsWith(">")) {
                        if (currChrom != "-1") {
                            // finished with this chromosome's sequence
                            println("fastaToNucSeq: finished chrom $currChrom")
                            chromNucSeqMap.put(currChrom, NucSeq(currSeq.toString()))
                        }
                        // reset chromosome name and sequence, begin processing next chrom
                        currChrom = line.substring(1).split(" ")[0]
                        currSeq = ByteArrayOutputStream()
                    } else {
                        currSeq.write(line.toByteArray())
                    }
                    line = br.readLine()
                }
                if (currSeq.size() > 0) {
                    println("fastaToNucSeq: finished chrom $currChrom")
                    chromNucSeqMap.put(currChrom, NucSeq(currSeq.toString()))
                }
            }
        } catch (exc:Exception) {
            myLogger.error("Error:  Exception in retrieveAgcData: ${exc.message}")
            throw IllegalStateException("Error:  Exception in retrieveAgcData: ${exc.message}")
        }

       // return agcOut
        return chromNucSeqMap
    }

}

