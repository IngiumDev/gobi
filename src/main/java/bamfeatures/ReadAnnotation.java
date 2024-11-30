package bamfeatures;

import gtf.structs.Gene;
import gtf.structs.Interval;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

import java.util.*;

public class ReadAnnotation {
    private TreeSet<Interval> firstRead;
    private TreeSet<Interval> secondRead;
    private String readID;
    private int splitCount;
    private int clippingSum;
    private int mismatchCount;
    private int geneCount;
    private int pcrIndex;
    private List<Gene> resultGenes;

    private boolean isFirstStrandNegative;
    private int alignmentStart;
    private int alignmentEnd;
    private boolean isReadFlipped;
    private int firstAlignmentStart;
    private int firstAlignmentEnd;
    private int secondAlignmentStart;
    private int secondAlignmentEnd;

    public ReadAnnotation() {
    }

    private static TreeSet<Interval> extractReadInterval(SAMRecord record) {
        TreeSet<Interval> intervals = new TreeSet<>();
        int lastEnd = -1;
        Interval currentRegion = null;
        for (AlignmentBlock block : record.getAlignmentBlocks()) {
            int start = block.getReferenceStart();
            int end = start + block.getLength() - 1;

            if (currentRegion == null || start > lastEnd + 1) {
                // Start a new region
                if (currentRegion != null) {
                    intervals.add(currentRegion);
                }
                currentRegion = new Interval(start, end);
            } else {
                // Merge with the existing region
                currentRegion.setEnd(end);
            }
            lastEnd = end;
        }

        // Add the last region if applicable
        if (currentRegion != null) {
            intervals.add(currentRegion);
        }
        return intervals;
    }

    public void extractReadAlignmentStartEnd(SAMRecord first, SAMRecord second) {

        // Reads can be flipped in the alignment
        isFirstStrandNegative = first.getReadNegativeStrandFlag();
        firstAlignmentStart = first.getAlignmentStart();
        firstAlignmentEnd = first.getAlignmentEnd();
        secondAlignmentStart = second.getAlignmentStart();
        secondAlignmentEnd = second.getAlignmentEnd();

        alignmentStart = Math.min(firstAlignmentStart, secondAlignmentStart);
        alignmentEnd = Math.max(firstAlignmentEnd, secondAlignmentEnd);

        // Determine if the read is flipped
        isReadFlipped = firstAlignmentStart >= secondAlignmentStart && (firstAlignmentStart != secondAlignmentStart || firstAlignmentEnd >= secondAlignmentEnd);

    }

    public void extractReadIntervals(SAMRecord first, SAMRecord second) {
        firstRead = extractReadInterval(first);
        secondRead = extractReadInterval(second);
    }


/* A read pair is
split inconsistent (annotate split-inconsistent: true) if there is at least one read base in one of the
reads that is within a split of the other read (skip all other annotations in these cases)
*/
public boolean areReadsConsistent() {
    // Check both directions since reads are interchangeable
    return checkConsistency(firstRead, secondRead) && checkConsistency(secondRead, firstRead);
}

    private boolean checkConsistency(TreeSet<Interval> readA, TreeSet<Interval> readB) {
        // Get the splits for readA
        List<Interval> splits = getSplits(readA);

        // Check every interval in readB against splits in readA
        for (Interval bInterval : readB) {
            for (Interval split : splits) {
                if (intervalsOverlap(bInterval, split)) {
                    return false; // Inconsistent
                }
            }
        }
        return true; // No inconsistency found
    }

    private List<Interval> getSplits(TreeSet<Interval> read) {
        List<Interval> splits = new ArrayList<>();
        Interval prev = null;

        for (Interval curr : read) {
            if (prev != null) {
                // Create a split interval between prev and curr
                int splitStart = prev.getEnd() + 1;
                int splitEnd = curr.getStart() - 1;
                if (splitStart <= splitEnd) { // Valid split
                    splits.add(new Interval(splitStart, splitEnd));
                }
            }
            prev = curr;
        }

        return splits;
    }

    private boolean intervalsOverlap(Interval a, Interval b) {
        return a.getStart() <= b.getEnd() && b.getStart() <= a.getEnd();
    }



}
