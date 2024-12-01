package gtf.treecollections;

import gtf.structs.Interval;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class StrandUnSpecificPCRIndex implements PCRIndexManager{
    Map<TreeSet<Interval>, Integer> pcrIndex;


    /**
     * Clears the PCR index and resets for the next chromosome
     */
    @Override
    public void nextChromosome() {
        pcrIndex = new HashMap<>();
    }

    /**
     * @param combinedRead                    the pair to be checked
     * @param isFirstStrandNegative
     * @return the PCR index of the pair
     */
    @Override
    public int getPCRIndex(TreeSet<Interval> combinedRead, boolean isFirstStrandNegative) {
        return pcrIndex.merge(combinedRead, 1, Integer::sum)-1;
    }

    @Override
    public void initializePCRIndex() {
        pcrIndex = new HashMap<>();
    }


}
