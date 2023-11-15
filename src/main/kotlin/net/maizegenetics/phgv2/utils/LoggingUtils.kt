package net.maizegenetics.phgv2.utils

import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import org.apache.logging.log4j.core.LoggerContext
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilderFactory

private val myOriginalOutputStream = System.out
private val myOriginalErrStream = System.err

/**
 * Sends all logging messages (including log4j) to standard out. The logging
 * level will be DEBUG.
 */
fun setupDebugLogging() {
    System.setOut(myOriginalOutputStream)
    System.setErr(myOriginalErrStream)
    sendDebugLog4jToStdout()
}

fun sendDebugLog4jToStdout() {

    val context = LogManager.getContext(false) as LoggerContext

    val builder = ConfigurationBuilderFactory.newConfigurationBuilder()

    val console = builder.newAppender("stdout", "Console")

    val standard = builder.newLayout("PatternLayout")
    standard.addAttribute("pattern", "[%t] %level %c{10} %d: %msg%n%throwable")

    console.add(standard)

    builder.add(console)

    val rootLogger = builder.newRootLogger(Level.DEBUG)
    rootLogger.add(builder.newAppenderRef("stdout"))
    builder.add(rootLogger)

    val configuration = builder.build()
    context.start(configuration)

}