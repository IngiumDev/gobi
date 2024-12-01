package gtf.treecollections;

import gtf.structs.Interval;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class StrandSpecificPCRIndex implements PCRIndexManager {
    Map<TreeSet<Interval>, Integer> plusStrandPCRIndex;
    Map<TreeSet<Interval>, Integer> minusStrandPCRIndex;

    /**
     * Clears the PCR index and resets for the next chromosome
     */
    @Override
    public void nextChromosome() {
        plusStrandPCRIndex = new HashMap<>();
        minusStrandPCRIndex = new HashMap<>();
    }
// PCR-index (format: pcrindex: ${pcrindex}): the number of read-pairs mapped exactly to
//the same genomic region as this read pair so far.

    /**
     * @param combinedRead          the pair to be checked
     * @param isFirstStrandNegative
     * @return the PCR index of the pair
     */
    @Override
    public int getPCRIndex(TreeSet<Interval> combinedRead, boolean isFirstStrandNegative) {
        if (isFirstStrandNegative) {
            return minusStrandPCRIndex.merge(combinedRead, 1, Integer::sum) - 1;
        } else {
            return plusStrandPCRIndex.merge(combinedRead, 1, Integer::sum) - 1;
        }
    }

    @Override
    public void initializePCRIndex() {
        // TODO: don't really need this as it is already initialized in next chromosome
        plusStrandPCRIndex = new HashMap<>();
        minusStrandPCRIndex = new HashMap<>();
    }
}
