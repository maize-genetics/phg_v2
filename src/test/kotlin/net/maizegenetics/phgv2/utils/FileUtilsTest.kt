package net.maizegenetics.phgv2.utils

import net.maizegenetics.phgv2.cli.TestExtension
import org.junit.jupiter.api.extension.ExtendWith
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import java.io.File
import java.io.FileNotFoundException
import kotlin.test.assertEquals

@ExtendWith(TestExtension::class)
class FileUtilsTest {

        @Test
        fun testGetBufferedWriter() {
            val tmpFile = File.createTempFile("test", ".txt")
            //use the path string because the method that takes a string calls the method that takes a file
            //  resulting in more complete testing coverage
            val tmpFilepath = tmpFile.absolutePath
            val bw = getBufferedWriter(tmpFilepath, false)
            bw.write(">Ref")
            bw.newLine()
            bw.write("ACGT")
            bw.newLine()
            bw.close()

            val br = getBufferedReader(tmpFilepath)
            val line = br.readLine()
            assertEquals(">Ref", line)
            val line2 = br.readLine()
            assertEquals("ACGT", line2)
        }

        @Test
        fun testGetBufferedWriterGzip() {
            val tmpFile = File.createTempFile("test", ".gz")
            val tmpFilepath = tmpFile.absolutePath
            val bw = getBufferedWriter(tmpFilepath, false)
            bw.write(">Ref")
            bw.newLine()
            bw.write("ACGT")
            bw.newLine()
            bw.close()

            val br = getBufferedReader(tmpFilepath)
            val line = br.readLine()
            assertEquals(">Ref", line)
            val line2 = br.readLine()
            assertEquals("ACGT", line2)
        }

        @Test
        fun testGetBufferedReaderURL() {
            val exception = assertThrows<FileNotFoundException> {
                getBufferedReader("https://s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa")
            }
            assertEquals(
                "https:/s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa (No such file or directory)",
                exception.message
            )
        }

        @Test
        fun testGetBufferedReaderURLGzip() {
            val exception = assertThrows<FileNotFoundException> {
                getBufferedReader("https://s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa.gz")
            }
            assertEquals(
                "https:/s3.amazonaws.com/maizegenetics/phg/phgV2Test/Ref.fa.gz (No such file or directory)",
                exception.message
            )
        }
}
