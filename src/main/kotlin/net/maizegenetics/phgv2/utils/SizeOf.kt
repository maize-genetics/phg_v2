package net.maizegenetics.phgv2.utils

//import net.maizegenetics.phgv2.utils.Sizeof.getMaxHeapSizeMB
//import net.maizegenetics.phgv2.utils.Sizeof.memoryUse
import java.text.NumberFormat

class SizeOf {
    companion object{
        fun printMemoryUse() {
            val format = NumberFormat.getInstance()
            try {
                println("-------------------------------")
                val current = memoryUse / 1048576L
                print("Current Heap Size: ")
                val currentStr = format.format(current)
                print(currentStr)
                println(" MB")
                print("Max Available Heap: ")
                System.out.print(getMaxHeapSizeMB())
                println(" MB")
                println("-------------------------------")
            } catch (e: Exception) {
                println("Problem getting heap size: " + e.message)
            }
        }
        /**
         * Returns max heap size in MB.
         *
         * @return max heap size
         */
        private fun getMaxHeapSizeMB(): Long {
            return Runtime.getRuntime().maxMemory() / 1048576L
        }

        private val memoryUse: Long
            get() {

                // Warm up all classes/methods we will use
                runGC()
                usedMemory()
                runGC()
                return usedMemory()
            }

        private fun runGC() {
            // It helps to call Runtime.gc()
            // using several method calls:
            for (i in 0..3) {
                _runGC()
            }
        }

        private fun _runGC() {
            var usedMem1 = usedMemory()
            var usedMem2 = Long.MAX_VALUE
            var i = 0
            while (usedMem1 < usedMem2 && i < 500) {
                s_runtime.runFinalization()
                s_runtime.gc()
                Thread.yield()
                usedMem2 = usedMem1
                usedMem1 = usedMemory()
                ++i
            }
        }

        private fun usedMemory(): Long {
            return s_runtime.totalMemory() - s_runtime.freeMemory()
        }

        private val s_runtime = Runtime.getRuntime()
    }


}
//object Sizeof {
//
//    @JvmStatic
//    fun main(args: Array<String>) {
//        printMemoryUse()
//    }
//
//    fun printMemoryUse() {
//        val format = NumberFormat.getInstance()
//        try {
//            println("-------------------------------")
//            val current = memoryUse / 1048576L
//            print("Current Heap Size: ")
//            val currentStr = format.format(current)
//            print(currentStr)
//            println(" MB")
//            print("Max Available Heap: ")
//            System.out.print(getMaxHeapSizeMB())
//            println(" MB")
//            println("-------------------------------")
//        } catch (e: Exception) {
//            println("Problem getting heap size: " + e.message)
//        }
//    }
//
//    /**
//     * Returns max heap size in MB.
//     *
//     * @return max heap size
//     */
//    private fun getMaxHeapSizeMB(): Long {
//        return Runtime.getRuntime().maxMemory() / 1048576L
//    }
//
//    private val memoryUse: Long
//        get() {
//
//            // Warm up all classes/methods we will use
//            runGC()
//            usedMemory()
//            runGC()
//            return usedMemory()
//        }
//
//    private fun runGC() {
//        // It helps to call Runtime.gc()
//        // using several method calls:
//        for (i in 0..3) {
//            _runGC()
//        }
//    }
//
//    private fun _runGC() {
//        var usedMem1 = usedMemory()
//        var usedMem2 = Long.MAX_VALUE
//        var i = 0
//        while (usedMem1 < usedMem2 && i < 500) {
//            s_runtime.runFinalization()
//            s_runtime.gc()
//            Thread.yield()
//            usedMem2 = usedMem1
//            usedMem1 = usedMemory()
//            ++i
//        }
//    }
//
//    private fun usedMemory(): Long {
//        return s_runtime.totalMemory() - s_runtime.freeMemory()
//    }
//
//    private val s_runtime = Runtime.getRuntime()
//}
