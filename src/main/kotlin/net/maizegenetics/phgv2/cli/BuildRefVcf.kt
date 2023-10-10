package net.maizegenetics.phgv2.cli

import biokotlin.genome.fastaToNucSeq
import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import com.google.common.collect.HashMultimap
import com.google.common.collect.Multimap
import com.google.common.collect.Range
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeaderLine
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.utils.*
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*
import java.util.logging.Logger


class BuildRefVcf : CliktCommand() {

    private val myLogger = Logger.getLogger("net.maizegenetics.phgv2.cli.BuildRefVcf")

    var myRefSequence: Map<String, NucSeq>? = null

    val bed by option(help = "BED formatted file of regions denoting reference haplotypes")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--bed must not be blank"
            }
        }

    val refFasta by option(help = "Reference FASTA file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--refFasta fasta must not be blank"
            }
        }

    val refName by option(help = "Line name for reference to be used in hvcf header")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference line name must not be blank"
            }
        }

    val outputDir by option("-o", "--output-dir", help = "Name for output VCF file Directory")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--output-dir/-o must not be blank"
            }
        }


    override fun run() {

        createRefHvcf(bed,refFasta,refName,outputDir)

    }

    fun createRefHvcf(ranges:String,refGenome:String,refName:String,outputDir:String) {

        myLogger.info("begin createRefHvcf,  refGenome=${refGenome}")


        // Verify the bed file is good.
        // If there are overlaps, throw an error and exit.
        // Overlaps are printed to stdout.
        val overlaps = verifyIntervalRanges(ranges)
        if (overlaps.isNotEmpty()) {
            // Overlaps not permitted.  User can fix via manually or via CreateValidIntervalsFilePlugin.  Throw error
            overlaps.stream().forEach { entry: String -> myLogger.severe("BuildRefVcf:  range Overlap entry: $entry") }
            throw IllegalArgumentException("BuildRefVcf: intervals bed file has overlapping positions. Overlapping reference ranges are not supported.  Please consolidate/remove overlaps.")
        }

        var groupAndPositionsMap: Multimap<String, Range<Position>> = HashMultimap.create()

        // read the ref fasta into the myRefSequence map
        myRefSequence = fastaToNucSeq(refGenome)

        // Create VCList: - this list contains VariantContext records for all the reference ranges
        val fullRefVCList = mutableListOf<VariantContext>()

        // Create HashSet to hold a list of ALT header lines for this gvcf
        // There will be 1 ALT header line created for every reference range in
        // this file
        val altHeaderLines: MutableSet<VCFHeaderLine> = HashSet()

        // Process the user interval ranges file
        try {
            bufferedReader(ranges).use { br ->
                var chrom = "-1"
                var prevChrom = "-1"
                var line: String? = null
                var chr: String? = null
                //var chrSeq = NucSeq("")
                var chrSeq: NucSeq? = NucSeq("")

                var chromAnchors = 0

                line = br.readLine()

                while (line != null) {
                    if (line.uppercase(Locale.getDefault()).contains("CHROMSTART")) {
                        line = br.readLine()
                        continue // skip header line
                    }
                    val tokens = line.split("\t")
                    // Line must contain at least 3 columns: chrom, chromStart, chromEnd
                    // Additional columns are ignored
                    if (tokens.size < 3) {
                        throw IllegalArgumentException("Error processing intervals file on line : ${line} . Must have values for columns chrom, chromStart, chromEnd and name")
                    }
                    chrom = tokens[0]

                    if (chrom != prevChrom) {
                        myLogger.info("Total intervals for chrom ${prevChrom}: ${chromAnchors}")
                        myLogger.info("Starting chrom $chrom")
                        chr = chrom
                        prevChrom = chrom
                        chromAnchors = 0
                        chrSeq = myRefSequence!![chr]
                    }

                    val anchorStart = tokens[1].toInt()  // NucSeq is 0 based, bed file is 0 based, no need to change
                    val anchorEnd = tokens[2].toInt()-1 // bed file is exclusive, decrement by 1 to make inclusive

                    chromAnchors++
                    // get bytes from reference, convert to string, add data to list
                    val intervalSeq = chrSeq!![anchorStart, anchorEnd].toString()
                    val intervalHash = getChecksumForString(intervalSeq, "Md5")
                    val intervalStart = Position(chrom, anchorStart)
                    val intervalEnd = Position(chrom, anchorEnd)
                    val intervalRange = Range.closed(intervalStart, intervalEnd)
                    val type = tokens[3]
                    groupAndPositionsMap.put(type, intervalRange)

                    val vc = createRefRangeVC(
                        myRefSequence!!,
                        refName,
                        intervalStart,
                        intervalEnd,
                        intervalStart,
                        intervalEnd
                    )
                    fullRefVCList.add(vc) // this is a list of ALL the VC records for all ranges - will become the hvcf file.

                    // Add vcf header lines here, doing somthing like this:
                    // headerLines.add(VCFAltHeaderLine("<ID=${intervalHash}, Description=\"${nodeDescription(node)}\">", VCFHeaderVersion.VCF4_2))
                    altHeaderLines.add(
                        VCFAltHeaderLine(
                            "<ID=${intervalHash}, Description=\"reference haplotype data\">,Number=6,Source=\"${refName}\",Contig=\"${chr}\",Start=\"${anchorStart}\",End=\"${anchorEnd}\",Checksum=\"Md5\",RefRange=\"${intervalHash}\">",
                            VCFHeaderVersion.VCF4_2
                        )
                    )
                    line = br.readLine()
                } // end while
                myLogger.info("Total intervals for chrom ${prevChrom} : ${chromAnchors}")
            } // end buffered reader

            // Load to an hvcf file, write files to user specified outputDir

            val hvcfFileName = "${refName}.hvcf"
            var localRefHVCFFile = outputDir + "/" + hvcfFileName

            //  This is in VariantUtils - it exports the gvcf file.
            // Include the VCF ALT Header lines created in the loop above
            exportVariantContext(refName, fullRefVCList, localRefHVCFFile, myRefSequence!!,altHeaderLines)
            //bgzip and csi index the file
            val bgzippedGVCFFileName = bgzipAndIndexGVCFfile(localRefHVCFFile)
            myLogger.info("${bgzippedGVCFFileName} created and stored to ${outputDir}")


        } catch (exc: Exception) {
            myLogger.severe("Error creating Ref HVCF file: ${exc.message}")
            throw IllegalStateException("Error creating Ref HVCF file: ${exc.message}")
        }

    } // end createRefRanges()

}