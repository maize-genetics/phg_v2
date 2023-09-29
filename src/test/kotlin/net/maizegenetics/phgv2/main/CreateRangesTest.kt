package net.maizegenetics.phgv2.main

import org.junit.jupiter.api.Test

class CreateRangesTest {
    @Test
    fun testCreateRangesCli() {
        val command = CreateRanges()

        command.run()

        assert(true)
    }
}