package net.maizegenetics.phgv2.main

import com.github.ajalt.clikt.testing.test
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals

class PhgTest {
    @Test
    fun testPhgCli() {
        val command = Phg()
        command.run()

        assert(true)
    }
}