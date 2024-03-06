package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.validate
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeaderLine
import net.maizegenetics.phgv2.utils.bgzipAndIndexGVCFfile
import net.maizegenetics.phgv2.utils.exportVariantContext
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * This class is used to create h.vcf files from existing PHG created g.vcf files
 * It only works on g.vcf files that were created by PHG and contain the assembly
 * and contig information in the header.
 * The db-path parameter is optional.  If it is not provided, the program will use
 * the current working directory.  It is expected this folder will contain the AGC
 * compressed genome data.
 *
 * Created h.vcf.gz and h.vcf.gz.csi files will be written to the same folder that contains the gvcf files.
 */
class Gvcf2Hvcf: CliktCommand(help = "Create  h.vcf files from existing PHG created g.vcf files")  {

    private val myLogger = LogManager.getLogger(Gvcf2Hvcf::class.java)

    val bed by option(help = "BED file with entries that define the haplotype boundaries")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--bed must not be blank"
            }
        }

    val referenceFile by option(help = "Path to local Reference FASTA file")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--reference-file must not be blank"
            }
        }

    val gvcfDir by option(help = "GVCF file directory. gvcf files should be bgzipped and csi indexed")
        .default("")
        .validate {
            require(it.isNotBlank()) {
                "--gvcf-dir must not be blank"
            }
        }

    val dbPath by option(help = "Folder name where TileDB datasets and AGC record is stored.  If not provided, the current working directory is used")
        .default("")

    override fun run() {
        // Check the dbPath and set it to the current working directory if it is not provided
        val dbPath = if (dbPath.isBlank()) {
            System.getProperty("user.dir")
        } else {
            dbPath
        }
        createASMHvcfs(dbPath, bed, referenceFile, gvcfDir)
    }

    // This is the entry point for this class.  It creates ranges from the bed file, stores the reference genome
    // sequence into a data structure, and then processes the gvcf files in the gvcfDir.  The hvcf files that
    // are created are written to the same folder that contains the gvcf files.
    fun createASMHvcfs(dbPath: String, bedFileName: String, referenceFileName: String, gvcfDirName: String) {
        //load the bed file into a data structure
        val ranges = CreateMafVcf().loadRanges(bedFileName)
        myLogger.info("CreateASMHvcfs: calling buildRefGenomeSeq")
        val refGenomeSequence = CreateMafVcf().buildRefGenomeSeq(referenceFileName)

        // walk the gvcf directory process files with g.vcf.gz extension
        File(gvcfDirName).walk().filter { !it.isHidden && !it.isDirectory }
            .filter { it.name.endsWith("g.vcf.gz")  || it.name.endsWith("g.vcf")  }.toList()
            .forEach {
                myLogger.info("CreateASMHvcfs: processing gvcf file: ${it.name}")
                // Create list of VariantContext from the gvcf file
                val vcfSampleAndVCs = createVCList(it)
                val asmHeaderLines = mutableMapOf<String, VCFHeaderLine>()
                val sampleName = vcfSampleAndVCs.first
                val variants = vcfSampleAndVCs.second
                //convert the GVCF records into hvcf records
                myLogger.info("createASMHvcfs: calling convertGVCFToHVCF for $sampleName")
                val hvcfVariants = CreateMafVcf().convertGVCFToHVCF(dbPath,sampleName, ranges, variants, refGenomeSequence, dbPath, asmHeaderLines)
                val asmHeaderSet = asmHeaderLines.values.toSet()
                //export the hvcfRecords
                myLogger.info("createASMHvcfs: calling exportVariantContext for $sampleName")
                var newFileName = it.name.replace(".g.vcf.gz", "")
                newFileName = newFileName.replace(".g.vcf", "")
                newFileName = "${newFileName}.h.vcf"

                exportVariantContext(sampleName, hvcfVariants, "${gvcfDirName}/${newFileName}",refGenomeSequence, asmHeaderSet)
                //bgzip the files
                bgzipAndIndexGVCFfile("${gvcfDirName}/${newFileName}")

            }
    }

    // This method takes a gvcfFile and returns a list of VariantContext records
    // It is assumed this is a single-sample gvcf file.
    fun createVCList(gvcfFile: File): Pair<String,List<VariantContext>> {
        // read the gvcf file and create a list of VariantContext records
        // for each line in the gvcf file, create a VariantContext record
        // and add it to the list
        // User should have both the bgzipped gvcf file PLUS an index for that file.
        // The index will be required later with API processing.  Here, we set requireIndex=true
        // so we can catch early cases where the index is not present.
        val gvcfFileReader = VCFFileReader(gvcfFile, false)

        // validate the gvcf file by checking if the INFO lines contain ID for
        // ASM_Chr and ASM_Start and ASM_End

        val infoHeaderLines = gvcfFileReader.fileHeader.getInfoHeaderLines()
        val asmChrFound = infoHeaderLines.any { it.id == "ASM_Chr" }
        val asmStartFound = infoHeaderLines.any { it.id == "ASM_Start" }
        val asmEndFound = infoHeaderLines.any { it.id == "ASM_End" }

        // This dies on the first non-PHG gvcf found.  It could make a list of all the non-PHG gvcfs
        // and just "continue" when a bad one is found.  But the exception makes it more apparent to the
        // user, who can then check remaining files that were not processed.
        if (!asmChrFound || !asmStartFound || !asmEndFound) {
            myLogger.error("The gvcf file ${gvcfFile} does not contain the required INFO lines: ASM_Chr, ASM_Start, ASM_End.\nThis file will cannot be processed.")
            throw IllegalStateException("The gvcf file does not contain the required INFO lines: ASM_Chr, ASM_Start, ASM_End")
        }

        var sampleNames = gvcfFileReader.fileHeader.sampleNamesInOrder
        val gvcfFileIterator = gvcfFileReader.iterator()

        val vcList = mutableListOf<VariantContext>()
        //Grab the first gvcf record
        var currentGVCFRecord:VariantContext? = gvcfFileIterator.next()

        while(currentGVCFRecord != null ) {
            vcList.add(currentGVCFRecord)
            currentGVCFRecord = gvcfFileIterator.next()
        }
        return Pair(sampleNames[0],vcList)

    }
}