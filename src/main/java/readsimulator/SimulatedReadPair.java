package readsimulator;

import gtf.structs.Exon;
import gtf.types.StrandDirection;

import java.util.SplittableRandom;
import java.util.TreeSet;

public class SimulatedReadPair extends Pair<SimulatedRead, SimulatedRead> {
    public SimulatedReadPair(SimulatedRead first, SimulatedRead second) {
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

    public SimulatedReadPair(String transcriptSequence, int fragmentStart, int fragmentLength, int readLength, String seqName, String geneID, String transcriptID, StrandDirection strandDirection, String reverseComplement) {
        super(new SimulatedRead(transcriptSequence, StrandDirection.FORWARD, fragmentStart, fragmentLength, readLength),
                new SimulatedRead(reverseComplement, StrandDirection.REVERSE, fragmentStart, fragmentLength, readLength));
        this.seqName = seqName;
        this.geneID = geneID;
        this.transcriptID = transcriptID;
        this.strandDirection = strandDirection;
    }

    public void mutateReadPairs(double mutationRate, SplittableRandom random, int mut1, int mut2) {
        this.getFirst().mutate(mutationRate, random, mut1);
        this.getSecond().mutate(mutationRate, random, mut2);
    }

    public void calculateGenomicPositions(TreeSet<Exon> exons) {
        this.getFirst().calculateGenomicPositions(exons, this.strandDirection);
        this.getSecond().calculateGenomicPositions(exons, this.strandDirection);
    }

}
