import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    kotlin("jvm") version "1.8.21"
    application
    id("org.jetbrains.kotlinx.kover") version "0.7.3"
}

group = "net.maizegenetics"
version = "2.0.0"

/*
This build script is need to use the early access
 */
buildscript {
    val kotlinVersion by extra("1.8.21")

    repositories {
        mavenCentral()
        gradlePluginPortal()
    }

    dependencies {
        classpath("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
        classpath(kotlin("serialization", version = kotlinVersion))
        classpath("org.jetbrains.dokka:dokka-gradle-plugin:1.6.21")
    }
}

repositories {
    mavenCentral()
    gradlePluginPortal()
    maven("https://maven.imagej.net/content/groups/public/")
}

dependencies {
    //testImplementation(kotlin("test"))
    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("org.biokotlin:biokotlin:0.08")

    implementation("com.github.ajalt.clikt:clikt:4.2.0")

    implementation("com.github.samtools:htsjdk:4.0.1")

    implementation("org.jetbrains.kotlin:kotlin-test:1.9.10")
    implementation("org.jetbrains.kotlin:kotlin-reflect:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-script-runtime:${kotlinVersion}")
    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.5.2")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.2.2")
    implementation("org.jetbrains.kotlinx:dataframe:0.8.0-rc-7")

    testImplementation("org.junit.jupiter:junit-jupiter:5.8.0")

    val kotestVersion = "5.6.2"
    listOf("runner-junit5", "assertions-core", "property", "framework-datatest").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }
}

tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
}

tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}

application {
    mainClass.set("main.PhgKt")

    // Set name of generated scripts in bin/
    applicationName = "phg"
}
