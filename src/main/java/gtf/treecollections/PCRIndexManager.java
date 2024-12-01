package gtf.treecollections;

import gtf.structs.Interval;

import java.util.TreeSet;

public interface PCRIndexManager {
    /**
     * Clears the PCR index and resets for the next chromosome
     */
    void nextChromosome();

    /**
     * @param combinedRead          the pair to be checked
     * @param isFirstStrandNegative
     * @return the PCR index of the pair
     */
    int getPCRIndex(TreeSet<Interval> combinedRead, boolean isFirstStrandNegative);

    void initializePCRIndex();
}
