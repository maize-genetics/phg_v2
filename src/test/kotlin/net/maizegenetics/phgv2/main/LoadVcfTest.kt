package net.maizegenetics.phgv2.main

import org.junit.jupiter.api.Test

class LoadVcfTest {
    @Test
    fun testLoadVcfCli() {
        val command = LoadVcf()

        command.run()

        assert(true)
    }
}