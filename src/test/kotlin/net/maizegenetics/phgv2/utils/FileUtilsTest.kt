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

            //test again with buffersize = 0
            val brBuffer0 = getBufferedReader(tmpFilepath, bufSize = 0)
            val line1Buffer0 = brBuffer0.readLine()
            assertEquals(">Ref", line1Buffer0)
            val line2Buffer0 = brBuffer0.readLine()
            assertEquals("ACGT", line2Buffer0)

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

}
