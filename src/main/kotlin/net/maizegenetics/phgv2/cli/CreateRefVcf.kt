package net.maizegenetics.phgv2.cli

import biokotlin.genome.fastaToNucSeq
import biokotlin.seq.NucSeq
import biokotlin.util.bufferedReader
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFHeaderLine
import htsjdk.variant.vcf.VCFHeaderVersion
import net.maizegenetics.phgv2.utils.*
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.util.*
import java.util.logging.Logger
import kotlin.io.path.nameWithoutExtension

/**
 * This class takes a reference genome and a bed file of intervals and creates a haplotype vcf for the reference genome.
 * The bed file and reference genome name are stored to the dbPath/reference folder
 * The hvcf file is stored to the dbPath/hvcf_files folder.
 * The hvcf file is loaded to tiledb.
 */
class CreateRefVcf : CliktCommand(help = "Create and load to tiledb a haplotype VCF for the reference genome") {

    private val myLogger = Logger.getLogger("net.maizegenetics.phgv2.cli.CreateRefVcf")

    var myRefSequence: Map<String, NucSeq>? = null

    // refurl is not required.  If present, it will result in a ##reference header
    // in the hvcf file.
    val referenceUrl by option(help = "URL where the reference FASTA file can be downloaded")
        .default("")

    val bed by option(help = "BED file with entries that define the haplotype boundaries")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--bed must not be blank"
            }
        }

    val referenceFile by option(help = "Path to local Reference FASTA file.")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }

    val referenceName by option(help = "Sample name for reference to be used in hvcf ")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-name must not be blank"
            }
        }

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    override fun run() {
        myLogger.info("begin CreateRefVcf: validate dpPath URI")
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }
        // Verify the dbPath is a valid tiledb URI
        val validDB = verifyURI(dbPath,"hvcf_dataset")

        createRefHvcf(bed,referenceFile,referenceName,referenceUrl,dbPath)
    }

    fun createRefHvcf(ranges:String,refGenome:String,refName:String,refUrl:String,outputDir:String) {

        myLogger.info("begin createRefHvcf,  refGenome=${refGenome}")


        // Verify the bed file is good.
        // If there are overlaps, throw an error and exit.
        // Overlaps are printed to stdout.
        val overlaps = verifyIntervalRanges(ranges)
        if (overlaps.isNotEmpty()) {
            // Overlaps not permitted.  User can fix via manually or via CreateValidIntervalsFilePlugin.  Throw error
            overlaps.forEach { entry: String -> myLogger.severe("BuildRefVcf:  range Overlap entry: $entry") }
            throw IllegalArgumentException("BuildRefVcf: intervals bed file has overlapping positions. Overlapping reference ranges are not supported.  Please consolidate/remove overlaps.")
        }

        // read the ref fasta into the myRefSequence map
        myRefSequence = fastaToNucSeq(refGenome)

        // Create VCList: - this list contains VariantContext records for all the reference ranges
        val fullRefVCList = mutableListOf<VariantContext>()

        // Create HashSet to hold a list of ALT header lines for this gvcf
        // There will be 1 ALT header line created for every reference range in
        // this file
        val altHeaderLines: MutableSet<VCFHeaderLine> = HashSet()

        // add ##reference line if a refURL was provided
        // The same is added for the assembly. When creating the reference hcvf,
        // we treat it like any other assembly vcf.  In this case, the reference and
        // assembly are the same.
        if (refUrl.isNotBlank()) {
            altHeaderLines.add(
                VCFHeaderLine("reference","${refUrl}")) // add the ##reference line
            altHeaderLines.add(
                VCFHeaderLine("assembly","${refUrl}")) // add the ##assembly line
        }

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
                        //chrSeq = myRefSequence!![chr]
                        if (myRefSequence!!.containsKey(chr)) {
                            chrSeq = myRefSequence!![chr]
                        } else {
                            throw IllegalStateException("Error processing intervals file on line : ${line} . Chromosome ${chr} not found in reference file")
                        }
                    }

                    val anchorStart = tokens[1].toInt() + 1  // NucSeq is 0 based, bed file is 0 based, but VCF is 1 based
                    val anchorEnd = tokens[2].toInt() // bed file is exclusive, but 0-based, so no need to change

                    chromAnchors++
                    // get bytes from reference, convert to string, add data to list
                    // And here it is tricky.  Bed file is 0-based, but we turned the start and end into 1-based
                    // to accommodate the VCF for later calls.  here we need to turn back to 0-based for NucSeq,
                    // which is inclusive/inclusive 0-based.  So the start is decremented but the end is not.
                    val intervalSeq = chrSeq!![anchorStart-1, anchorEnd-1].toString()
                    val intervalHash = getChecksumForString(intervalSeq, "Md5")
                    val intervalStart = Position(chrom, anchorStart)
                    val intervalEnd = Position(chrom, anchorEnd)

                    // create the calls Pair.  This is a Pair of the reference allele and the interval hash
                    // The interval hash becomes the alt allele in the hvcf file
                    val refAllele = chrSeq!![anchorStart-1].toString() // subtract 1 to convert to 0-based

                    val calls = Pair(refAllele,intervalHash)
                    // this prevents them from being written to the hvcf file
                    // The intervalStart/intervalEnd values passed here are 1-based
                    val vc = createHVCFRecord(
                        refName,
                        intervalStart,
                        intervalEnd,
                        calls
                    )

                    fullRefVCList.add(vc) // this is a list of ALL the VC records for all ranges - will become the hvcf file.

                    // Add vcf header lines here, doing somthing like this:
                    // headerLines.add(VCFAltHeaderLine("<ID=${intervalHash}, Description=\"${nodeDescription(node)}\">", VCFHeaderVersion.VCF4_2))
                    altHeaderLines.add(
                        VCFAltHeaderLine(
                              "<ID=${intervalHash}, Description=\"haplotype data for line: ${refName}\"," +
                                       "Source=\"${refGenome}\",SampleName=\"${refName}\",Regions=\"${chr}:${anchorStart}-${anchorEnd}\"," +
                                        "Checksum=\"Md5\",RefRange=\"${intervalHash}\">",
                        VCFHeaderVersion.VCF4_2
                    )

                    )
                    line = br.readLine()
                } // end while
                myLogger.info("Total intervals for chrom ${prevChrom} : ${chromAnchors}")
            } // end buffered reader

            // Load to an hvcf file, write files to user specified outputDir

            val hvcfFileName = "${refName}.h.vcf"
            val tiledbHvcfDir = "${dbPath}/hvcf_files"

            // Make a folder for the hvcf files
            // ALl hvcf files loaded to this tiledb dataset will be saved in this subfolder
            Files.createDirectories(Paths.get(tiledbHvcfDir)) // skips creating folders that already exist
            var localRefHVCFFile = tiledbHvcfDir + "/" + hvcfFileName

            //  This is in VariantUtils - it exports the gvcf file.
            // Include the VCF ALT Header lines created in the loop above
            exportVariantContext(refName, fullRefVCList, localRefHVCFFile, myRefSequence!!,altHeaderLines)
            //bgzip and csi index the file
            val bgzippedGVCFFileName = bgzipAndIndexGVCFfile(localRefHVCFFile)
            myLogger.info("${bgzippedGVCFFileName} created and stored to ${tiledbHvcfDir}")

            // Create the tiledb datasets if they don't exist
            File(dbPath).mkdirs()
            Initdb().createDataSets(dbPath)

            // Load the Ref hvcf file to tiledb
            val loadVcf = LoadVcf()
            loadVcf.loadVcfFiles(tiledbHvcfDir,dbPath) // default the threads to "1" for the ref.

            // Store the bed file and the reference fasta file to the dbPath/reference folder
            // Create the dbPath/reference folder if it doesn't exist
            val tiledbRefDir = "${dbPath}/reference"
            Files.createDirectories(Paths.get(tiledbRefDir)) // skips folders that already exist
            val fileName = Paths.get(ranges).nameWithoutExtension
            val localBedFile = "${tiledbRefDir}/${fileName}.bed"

            // The ref fasta file will be stored as <refName>.fa.  This allows methods
            // in other parts of the pipeline to know the sampleName of the reference genome as it
            // is stored in tileDB and AGC.
            val localRefFile = "${tiledbRefDir}/${refName}.fa"
            Files.copy(Paths.get(ranges), Paths.get(localBedFile))
            Files.copy(Paths.get(refGenome), Paths.get(localRefFile))

        } catch (exc: Exception) {
            myLogger.severe("Error creating Ref HVCF file: ${exc.message}")
            throw IllegalStateException("Error creating Ref HVCF file: ${exc.message}")
        }

    } // end createRefRanges()

}