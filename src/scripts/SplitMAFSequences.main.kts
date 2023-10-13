import java.io.File

/**
 * This script splits the sequences in a MAF file into smaller blocks.
 * The input argument in the original MAF file. The output file will
 * be named <input filename>.split
 */

val filename = args[0]
println("input: $filename")

val output = "$filename.split"
println("output: $output")

val maxSequenceLength = 300

File(filename).bufferedReader().use { reader ->

    File(output).bufferedWriter().use { writer ->

        var line = reader.readLine()
        while (line != null) {

            when {

                line.isBlank() -> {
                    // skip
                    line = reader.readLine()
                }

                line.startsWith("##") -> {
                    writer.write(line)
                    writer.write("\n")
                    println("header: $line")
                    line = reader.readLine()
                }

                line.startsWith("a") -> {

                    val alignBlock = line

                    line = reader.readLine()
                    val sequences = mutableListOf<String>()
                    while (line != null && line.startsWith("s")) {
                        sequences.add(line)
                        line = reader.readLine()
                    }

                    val numSequences = sequences.size
                    val src = mutableListOf<String>()
                    val start = mutableListOf<String>()
                    val size = mutableListOf<String>()
                    val strand = mutableListOf<String>()
                    val srcSize = mutableListOf<String>()
                    val sequence = mutableListOf<String>()
                    sequences.forEach { line ->
                        val tokens = line.split("\t")
                        src.add(tokens[1])
                        start.add(tokens[2])
                        size.add(tokens[3])
                        strand.add(tokens[4])
                        srcSize.add(tokens[5])
                        sequence.add(tokens[6])
                    }

                    val subseqs = mutableListOf<List<String>>()
                    sequence.forEach {
                        subseqs.add(it.windowed(maxSequenceLength, maxSequenceLength, true))
                    }

                    val currentStart = mutableListOf<Int>()
                    for (i in 0 until numSequences) {
                        println("start: ${start[i]}")
                        currentStart.add(start[i].trim().toInt())
                    }

                    val numSubSeqs = subseqs[0].size
                    for (i in 0 until numSubSeqs) {
                        writer.write(alignBlock)
                        writer.write("\n")
                        for (j in 0 until numSequences) {
                            val subseqLenNoDashes = subseqs[j][i].filter { it != '-' }.length
                            writer.write("s\t")
                            writer.write(src[j])
                            writer.write("\t")
                            writer.write(currentStart[j].toString())
                            writer.write("\t")
                            writer.write(subseqLenNoDashes.toString())
                            writer.write("\t")
                            writer.write(strand[j])
                            writer.write("\t")
                            writer.write(srcSize[j])
                            writer.write("\t")
                            writer.write(subseqs[j][i])
                            writer.write("\n")

                            currentStart[j] += subseqLenNoDashes
                        }

                        writer.write("\n")

                    }

                }

                else -> {
                    throw Exception("Unexpected line: $line")
                }
            }

        }

    }

}