import java.io.File

val inputFile =
    "/Users/terry/hackathon/20250127-ortholog-contrast/zero-shot-genome-browser/1001_genomes_variants_zero_shot_first1000.tsv"
val outputFile =
    "/Users/terry/hackathon/20250127-ortholog-contrast/zero-shot-genome-browser/1001_genomes_variants_zero_shot_first1000.wig"

File(outputFile).bufferedWriter().use { writer ->

    writer.write("track type=wiggle_0 name=\"My Data\" description=\"My Data Description\"\n")
    var currentChromosome: String? = null

    File(inputFile).forEachLine { line ->

        // Skip header line
        if (line.startsWith("chr\tpos\t")) return@forEachLine

        // Split line into components
        val columns = line.split("\t")
        val chromosome = columns[0]
        val position = columns[1]
        val value = columns[2]

        // If the chromosome changes, write a new variableStep header
        if (chromosome != currentChromosome) {
            writer.write("variableStep chrom=$chromosome\n")
            currentChromosome = chromosome
        }

        // Write position and value to the output file
        writer.write("$position $value\n")
    }

}
