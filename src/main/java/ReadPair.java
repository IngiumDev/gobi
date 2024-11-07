import gtf.structs.Interval;
import gtf.types.StrandDirection;

import java.util.Random;
import java.util.TreeSet;

public class ReadPair extends Pair<Read, Read> {
    public ReadPair(Read first, Read second) {
        super(first, second);
    }


    public ReadPair(String transcriptSequence, int fragmentStart, int fragmentLength, int readLength) {
        super(new Read(transcriptSequence, StrandDirection.FORWARD, fragmentStart, fragmentLength, readLength),
                new Read(transcriptSequence, StrandDirection.REVERSE, fragmentStart, fragmentLength, readLength));
    }

    public void mutateReadPairs(double mutationRate, Random random) {
        this.getFirst().mutate(mutationRate, random);
        this.getSecond().mutate(mutationRate, random);
    }

    public void calculateGenomicPositions(TreeSet<Interval> exonPositions) {
        this.getFirst().calculateGenomicPositions(exonPositions);
        this.getSecond().calculateGenomicPositions(exonPositions);
    }

}
