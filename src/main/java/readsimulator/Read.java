package readsimulator;

import gtf.structs.Interval;
import gtf.types.StrandDirection;
import parsers.GenomeSequenceExtractor;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import static runners.ReadSimulatorRunner.getCoveredRegion;

public class Read {
    public final static char[] nucleotides = {'A', 'C', 'G', 'T'};
    private final Interval transcriptCoordinates;
    private final StrandDirection strandDirection;
    private String seq;
    private List<Integer> mutatedPositions;
    private TreeSet<Interval> chromosomalCoordinates;

    public Read(String transcriptSeq, StrandDirection direction, int fragmentStart, int fragmentLength, int readLength) {
        if (direction == StrandDirection.FORWARD) {
            this.seq = transcriptSeq.substring(fragmentStart, fragmentStart + readLength);
            this.transcriptCoordinates = new Interval(fragmentStart, fragmentStart + readLength - 1);
            this.strandDirection = direction;
        } else if (direction == StrandDirection.REVERSE) {
            this.seq = GenomeSequenceExtractor.reverseComplement(transcriptSeq.substring(fragmentStart + fragmentLength - readLength, fragmentStart + fragmentLength));
            this.transcriptCoordinates = new Interval(fragmentStart + fragmentLength - readLength, fragmentStart + fragmentLength - 1);
            this.strandDirection = direction;
        } else {
            throw new IllegalArgumentException("StrandDirection not given");
        }
    }

    public static void main(String[] args) {
        String seq = "ATTTTTA";
        int fragmentStart = 0;
        int fragmentLength = 7;
        int readLength = 4;
        Read read1 = new Read(seq, StrandDirection.FORWARD, fragmentStart, fragmentLength, readLength);
        Read read2 = new Read(seq, StrandDirection.REVERSE, fragmentStart, fragmentLength, readLength);
        System.out.println();
    }

    public Interval getTranscriptCoordinates() {
        return transcriptCoordinates;
    }

    public StrandDirection getStrandDirection() {
        return strandDirection;
    }

    public String getSeq() {
        return seq;
    }

    public List<Integer> getMutatedPositions() {
        return mutatedPositions;
    }

    public TreeSet<Interval> getChromosomalCoordinates() {
        return chromosomalCoordinates;
    }

    // TODO: change so that mutation's have to happen 33% chance for each nucleotide
    public void mutate(double mutationRate, Random random) {
        StringBuilder mutatedSeq = new StringBuilder(seq);
        mutatedPositions = new ArrayList<>();
        for (int i = 0; i < seq.length(); i++) {
            if (random.nextDouble() < mutationRate) {
                mutatedSeq.setCharAt(i, getRandomNucleotide(mutatedSeq.charAt(i), random));
                mutatedPositions.add(i);
            }
        }
        this.seq = mutatedSeq.toString();
    }

    public char getRandomNucleotide(char currentNucleotide, Random random) {
        char newNucleotide;
        do {
            newNucleotide = nucleotides[random.nextInt(nucleotides.length)];
        } while (newNucleotide == currentNucleotide);
        return newNucleotide;
    }

    public void calculateGenomicPositions(TreeSet<Interval> exonPositions) {
        chromosomalCoordinates = getCoveredRegion(exonPositions, transcriptCoordinates);
    }
}
