package readsimulator;

import gtf.structs.Interval;
import gtf.types.StrandDirection;
import org.apache.commons.rng.UniformRandomProvider;

import java.util.Random;
import java.util.TreeSet;

public class ReadPair extends Pair<Read, Read> {
    public ReadPair(Read first, Read second) {
        super(first, second);
    }
    private String seqName;
    private String geneID;
    private String transcriptID;
    private StrandDirection strandDirection;

    public String getSeqName() {
        return seqName;
    }

    public String getGeneID() {
        return geneID;
    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public ReadPair(String transcriptSequence, int fragmentStart, int fragmentLength, int readLength, String seqName, String geneID, String transcriptID, StrandDirection strandDirection) {
        super(new Read(transcriptSequence, StrandDirection.FORWARD, fragmentStart, fragmentLength, readLength),
                new Read(transcriptSequence, StrandDirection.REVERSE, fragmentStart, fragmentLength, readLength));
        this.seqName = seqName;
        this.geneID = geneID;
        this.transcriptID = transcriptID;
        this.strandDirection = strandDirection;
    }

    public void mutateReadPairs(double mutationRate, Random random, UniformRandomProvider rng) {
        this.getFirst().mutate(mutationRate, random,rng);
        this.getSecond().mutate(mutationRate, random,rng);
    }

    public void calculateGenomicPositions(TreeSet<Interval> exonPositions) {
        this.getFirst().calculateGenomicPositions(exonPositions, this.strandDirection);
        this.getSecond().calculateGenomicPositions(exonPositions, this.strandDirection);
    }

}
