import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

val ktorVersion = "2.3.7"

plugins {
    kotlin("jvm") version "1.9.24"
    application
    id("org.jetbrains.kotlinx.kover") version "0.7.3"
    kotlin("plugin.serialization") version "1.6.21"
}

group = "net.maizegenetics"

/*
This build script is need to use the early access
 */
buildscript {
    val kotlinVersion by extra("1.9.24")

    repositories {
        mavenCentral()
        gradlePluginPortal()
    }

    dependencies {
        classpath("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
        classpath(kotlin("serialization", version = kotlinVersion))
        classpath("org.jetbrains.dokka:dokka-gradle-plugin:1.9.20")
    }
}

repositories {
    mavenCentral()
    gradlePluginPortal()
    maven("https://maven.imagej.net/content/groups/public/")
}

dependencies {

    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("org.jetbrains.kotlinx:dataframe:0.13.1")
    implementation("org.jetbrains.lets-plot:lets-plot-kotlin:4.7.0")

    implementation("org.jetbrains.lets-plot:lets-plot-kotlin-jvm:4.7.0")
    implementation("org.jetbrains.lets-plot:lets-plot-image-export:4.3.0")
    implementation("org.jetbrains.lets-plot:lets-plot-batik:4.3.0")
    implementation("org.biokotlin:biokotlin:0.26")
    implementation("com.github.ajalt.clikt:clikt:4.2.0")

    implementation("com.github.samtools:htsjdk:4.0.1")

    implementation("it.unimi.dsi:fastutil:8.5.12")
    implementation("net.openhft:zero-allocation-hashing:0.16")

    implementation("org.jetbrains.kotlin:kotlin-test:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-reflect:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-script-runtime:${kotlinVersion}")
    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.8.1")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.2.2")

    implementation("org.apache.commons:commons-math3:3.6.1")
    implementation("org.apache.logging.log4j:log4j-core:2.20.0")
    implementation("org.apache.logging.log4j:log4j-api:2.20.0")

    implementation("com.google.guava:guava:33.1.0-jre")

    // Use these jar file when compiling for Linux
    // Keep the Mac Intel and Mac ARM jar inclusions commented out
    implementation(files("repo/tiledb-vcf-java-0.37.0-1-g03553439.jar")) // added 2024-01-08 based on TileDB VCF 0.37.0 build
    // TILedb universal Java API
    implementation("io.tiledb:tiledb-java:0.28.1") // Replace with the latest version
    //implementation(files("repo/tiledb-java-0.19.6-SNAPSHOT.jar"))

    // Use these jar files when compiling for Mac with Intel chip
    // Keep the Linux and Mac ARM jar inclusions commented out
//    implementation(files("repo/MacIntel_tiledb-vcf-java-0.28.0.jar"))
//    implementation(files("repo/MacIntel_tiledb-java-0.21.1-SNAPSHOT.jar"))

    // Use these jar files when compiling for Mac with ARM chip
    // Keep the Linux and Mac Intel jar inclusions commented out
//    implementation(files("repo/MacARM_tiledb-vcf-java-0.28.0.jar"))
//    implementation(files("repo/MacARM_tiledb-java-0.21.1-SNAPSHOT.jar"))

    implementation("io.ktor:ktor-server-netty:$ktorVersion")
    implementation("io.ktor:ktor-server-content-negotiation:$ktorVersion")
    implementation("io.ktor:ktor-serialization-kotlinx-json:$ktorVersion")
    implementation("io.ktor:ktor-server-websockets:$ktorVersion")
    implementation("io.ktor:ktor-server-default-headers:$ktorVersion")
    implementation("io.ktor:ktor-server-call-logging:$ktorVersion")
    implementation("io.ktor:ktor-server-content-negotiation:$ktorVersion")
    implementation("io.ktor:ktor-serialization-jackson:$ktorVersion")
    implementation("io.ktor:ktor-client-core:$ktorVersion")
    implementation("io.ktor:ktor-client-cio:$ktorVersion")
    implementation("io.ktor:ktor-client-content-negotiation:$ktorVersion")
    implementation("io.ktor:ktor-client-websockets:$ktorVersion")

    testImplementation("org.junit.jupiter:junit-jupiter:5.8.0")


    // for testing ktor
    testImplementation("io.ktor:ktor-server-test-host:$ktorVersion")
    testImplementation("io.ktor:ktor-client-mock:$ktorVersion")
    testImplementation("org.jetbrains.kotlin:kotlin-test:$kotlinVersion")


    val kotestVersion = "5.6.2"
    listOf("runner-junit5", "assertions-core", "property", "framework-datatest").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }
}

// include versions.properties file in jar file
tasks.jar {
    from(sourceSets.main.get().output)
    from(projectDir) {
        include("version.properties")
    }
    exclude("application.conf")
}

tasks.distTar {
    from("${projectDir}/src/main/resources/application.conf") {
        into("phg/resources/main/")
    }
}

tasks.startScripts {
    classpath = classpath?.plus(files("resources"))
    doLast {
        val windowsScriptFile = windowsScript
        val unixScriptFile = unixScript

        windowsScriptFile.writeText(
            windowsScriptFile.readText().replace("%APP_HOME%\\lib\\resources", "%APP_HOME%\\resources\\main")
        )
        unixScriptFile.writeText(
            unixScriptFile.readText().replace("\$APP_HOME/lib/resources", "\$APP_HOME/resources/main")
        )
    }
}

koverReport {
    verify {
        rule {
            "Minimal line coverage rate as a percentage"
            bound {
                minValue = 70
            }
        }
    }
}


tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
//    Uncomment this if you want to increase the heap size for testing
//    maxHeapSize = "10G"
}

tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "17"
}

application {
    mainClass.set("net.maizegenetics.phgv2.cli.PhgKt")

    // Set name of generated scripts in bin/
    applicationName = "phg"
}
