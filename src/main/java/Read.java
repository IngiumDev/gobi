import gtf.structs.Interval;
import gtf.types.StrandDirection;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

public class Read {
    public final static char[] nucleotides = {'A', 'C', 'G', 'T'};
    private String seq;
    private final Interval transcriptCoordinates;
    private final StrandDirection strandDirection;
    private List<Integer> mutatedPositions;
    private TreeSet<Interval> chromosomalCoordinates;

    public Read(String seq, Interval transcriptCoordinates, StrandDirection strandDirection) {
        this.seq = seq;
        this.transcriptCoordinates = transcriptCoordinates;
        this.strandDirection = strandDirection;
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
}
