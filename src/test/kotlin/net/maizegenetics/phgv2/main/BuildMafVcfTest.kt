package net.maizegenetics.phgv2.main

import org.junit.jupiter.api.Test

class BuildMafVcfTest {
    @Test
    fun testBuildMafVcfCli() {
        val command = BuildMafVcf()

        command.run()

        assert(true)
    }
}