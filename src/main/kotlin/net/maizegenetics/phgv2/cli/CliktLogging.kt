package net.maizegenetics.phgv2.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.OptionWithValues
import org.apache.logging.log4j.LogManager

private val myLogger = LogManager.getLogger("net.maizegenetics.phgv2.cli.CliktLogging")


/**
 * Log the command name and all parameters and their values.
 */
fun logCommand(command: CliktCommand) {

    command.currentContext.originalArgv.joinToString(" ").let {
        myLogger.info("PHGv2 Command: $it")
    }

    command.currentContext.command.registeredOptions().forEach {
        try {
            val option = it as OptionWithValues<*, *, *>
            println("    ${option.names} = ${option.value}")
        } catch (e: Exception) {
            // do nothing
        }
    }

}
fun headerCommand(command: CliktCommand) : String {
    return "PHGv2 Command: ${command.currentContext.originalArgv.joinToString(" ")}"
}
