package net.maizegenetics.phgv2.main

import org.junit.jupiter.api.Test

class InitdbTest {
    @Test
    fun testInitdbCli() {
        val command = Initdb()

        command.run()

        assert(true)
    }
}