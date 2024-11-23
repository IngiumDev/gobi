package readsimulator;

import gtf.structs.Exon;
import gtf.structs.Interval;
import gtf.types.StrandDirection;
import parsers.GenomeSequenceExtractor;

import java.util.HashSet;
import java.util.Set;
import java.util.SplittableRandom;
import java.util.TreeSet;

import static runners.ReadSimulatorRunner.getCoveredRegion;

public class Read {
    public final static char[] NUCLEOTIDES
            = {'A', 'C', 'G', 'T'};
    private final Interval transcriptCoordinates;
    private final StrandDirection strandDirection;
    private final int transcriptSeqLength;
    private String seq;
    private Set<Integer> mutatedPositions;
    private TreeSet<Interval> chromosomalCoordinates;

    public Read(String transcriptSeq, StrandDirection direction, int fragmentStart, int fragmentLength, int readLength) {


        if (direction == StrandDirection.FORWARD) {
            this.seq = transcriptSeq.substring(fragmentStart, fragmentStart + readLength);
            this.transcriptCoordinates = new Interval(fragmentStart, fragmentStart + readLength - 1);
            this.strandDirection = direction;
            this.transcriptSeqLength = transcriptSeq.length();
        } else if (direction == StrandDirection.REVERSE) {
            this.seq = GenomeSequenceExtractor.reverseComplement(transcriptSeq.substring(fragmentStart + fragmentLength - readLength, fragmentStart + fragmentLength));
            this.transcriptCoordinates = new Interval(fragmentStart + fragmentLength - readLength, fragmentStart + fragmentLength - 1);
            this.strandDirection = direction;
            this.transcriptSeqLength = transcriptSeq.length();
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

    public Set<Integer> getMutatedPositions() {
        return mutatedPositions;
    }

    public TreeSet<Interval> getChromosomalCoordinates() {
        return chromosomalCoordinates;
    }

    public void mutate(double mutationRate, SplittableRandom random, int numMutations) {
        mutatedPositions = new HashSet<>();
        if (numMutations == 0) {
            return;
        }
        StringBuilder mutatedSeq = new StringBuilder(seq);
        int length = seq.length();
        int mutated = 0;
        while (mutated < numMutations) {

            int position = random.nextInt(length);
            if (mutatedPositions.contains(position)) {
                continue;
            }

            mutatedPositions.add(position);
            mutatedSeq.setCharAt(position, getRandomNucleotide(seq.charAt(position), random));
            mutated++;
        }
        seq = mutatedSeq.toString();
    }

    private int sampleMutations(SplittableRandom random, int readLength, double mutationRate) {
        int numMutations = 0;
        for (int i = 0; i < readLength; i++) {
            if (random.nextDouble() < mutationRate) {
                numMutations++;
            }
        }
        return numMutations;
    }


    public char getRandomNucleotide(char currentNucleotide, SplittableRandom random) {
        char newNucleotide;
        do {
            newNucleotide = NUCLEOTIDES[random.nextInt(NUCLEOTIDES.length)];
        } while (newNucleotide == currentNucleotide);
        return newNucleotide;
    }

    public void calculateGenomicPositions(TreeSet<Exon> exons, StrandDirection direction) {
        if (direction == StrandDirection.FORWARD) {

            chromosomalCoordinates = getCoveredRegion(exons, transcriptCoordinates, direction);

        } else {
            Interval reversedCoordinates = new Interval(transcriptSeqLength - transcriptCoordinates.getEnd() - 1, transcriptSeqLength - transcriptCoordinates.getStart() - 1);
            chromosomalCoordinates = getCoveredRegion(exons, reversedCoordinates, direction);
        }
    }
}
