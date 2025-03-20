package net.maizegenetics.phgv2.brapi.service

import net.maizegenetics.phgv2.brapi.model.Sample
import net.maizegenetics.phgv2.brapi.utilities.BrAPIConfig
import net.maizegenetics.phgv2.utils.inputStreamProcessing
import org.apache.logging.log4j.LogManager
import java.io.BufferedInputStream

object SamplesService {

    private val myLogger = LogManager.getLogger(SamplesService::class.java)

    // Cached map of all taxa. Key is genoid mapped to Sample object
    private val taxa: Map<String, Sample> by lazy {
        taxaMap("${BrAPIConfig.tiledbURI}/hvcf_dataset")
    }

    private val taxaByName: Map<String, Sample> by lazy {
        taxa.values.associateBy { it.sampleName }
    }

    fun allTaxaNames(): List<Sample> {
        return taxa.values.toList()
    }

    // This function runs tiledbvcf list --uri <uri> and returns a map of sample names
    private fun taxaMap(uri: String): Map<String, Sample> {
        try {
            var builder = ProcessBuilder("conda", "run", "-n", "phgv2-tiledb", "tiledbvcf", "list", "--uri", uri)
                .start()
            val sampleListOut = BufferedInputStream(builder.inputStream, 5000000)
            val temp = inputStreamProcessing(sampleListOut)
            return temp.mapIndexed { index, sample -> Sample(sample, index.toString()) }.associateBy { it.sampleDbId }

        } catch (exc: Exception) {
            myLogger.error("Error reading tiledb list output on ${uri}.")
            throw IllegalArgumentException("Error running ProcessBuilder to list samples from ${uri}: $exc")
        }
    }

    fun taxa(id: String) = taxa[id]

}