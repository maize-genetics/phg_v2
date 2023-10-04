package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.context
import com.github.ajalt.clikt.core.subcommands
import com.github.ajalt.clikt.parameters.options.versionOption
import com.github.ajalt.mordant.rendering.AnsiLevel
import com.github.ajalt.mordant.terminal.Terminal

class Phg : CliktCommand() {

    // Need an automated way to get version from build.gradle.kts
    private val version = "2.0.0"

    init {
        versionOption(version)
    }

    override fun run() = Unit

}

fun main(args: Array<String>) = Phg()
    .subcommands(Initdb(), CreateRanges(), AlignAssemblies(), BuildRefVcf(), BuildMafVcf(), LoadVcf())
    .context {
        terminal = Terminal(ansiLevel = AnsiLevel.TRUECOLOR, interactive = true)
    }.main(args)