package gtf.structs;

public class GTFTimer {
    static long gtfParseTime = 0;
    static long intronProcessTime = 0;
    static long exonProcessTime = 0;
    static long outputTime = 0;
    static long totalTime = 0;

    public static long getIntronProcessTime() {
        return intronProcessTime;
    }

    public static void setIntronProcessTime(long intronProcessTime) {
        GTFTimer.intronProcessTime = intronProcessTime;
    }

    public static long getGtfParseTime() {
        return gtfParseTime;
    }

    public static void setGtfParseTime(long gtfParseTime) {
        GTFTimer.gtfParseTime = gtfParseTime;
    }

    public static long getExonProcessTime() {
        return exonProcessTime;
    }

    public static void setExonProcessTime(long exonProcessTime) {
        GTFTimer.exonProcessTime = exonProcessTime;
    }

    public static long getOutputTime() {
        return outputTime;
    }

    public static void setOutputTime(long outputTime) {
        GTFTimer.outputTime = outputTime;
    }

    public static long getTotalTime() {
        return totalTime;
    }

    public static void setTotalTime(long totalTime) {
        GTFTimer.totalTime = totalTime;
    }
}
