package readsimulator;

import gtf.structs.Exon;
import gtf.structs.Interval;
import gtf.types.StrandDirection;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.statistics.distribution.BinomialDistribution;
import parsers.GenomeSequenceExtractor;

import java.util.*;

import static runners.ReadSimulatorRunner.getCoveredRegion;

public class Read {
    public final static char[] NUCLEOTIDES
            = {'A', 'C', 'G', 'T'};
    private final Interval transcriptCoordinates;
    private final StrandDirection strandDirection;
    private final int transcriptSeqLength;
    private String seq;
    private List<Integer> mutatedPositions;
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

    public List<Integer> getMutatedPositions() {
        return mutatedPositions;
    }

    public TreeSet<Interval> getChromosomalCoordinates() {
        return chromosomalCoordinates;
    }

    // TODO: Optimize
    public void mutate(double mutationRate, Random random, UniformRandomProvider rng) {
        BinomialDistribution binomialDistribution = BinomialDistribution.of(seq.length(), mutationRate);
        int numMutations = binomialDistribution.createSampler(rng).sample();
        mutatedPositions = new ArrayList<>(numMutations);

        if (numMutations == 0) {
            return;
        }

        StringBuilder mutatedSeq = new StringBuilder(seq);
        Set<Integer> positionsToMutate = new HashSet<>();
        while (positionsToMutate.size() < numMutations) {
            positionsToMutate.add(random.nextInt(seq.length()));
        }

        // Apply mutations at the unique positions
        mutatedPositions = new ArrayList<>(positionsToMutate); // Store mutation positions
        for (int position : positionsToMutate) {
            mutatedSeq.setCharAt(position, getRandomNucleotide(mutatedSeq.charAt(position), random));
        }

        // Convert the mutated char array back to a string
        this.seq = mutatedSeq.toString();

    }


    public char getRandomNucleotide(char currentNucleotide, Random random) {
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
