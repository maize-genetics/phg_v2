package net.maizegenetics.phgv2.main

import org.junit.jupiter.api.Test

class BuildRefVcfTest {
    @Test
    fun testBuildRefVcfCli() {
        val command = BuildRefVcf()

        command.run()

        assert(true)
    }
}