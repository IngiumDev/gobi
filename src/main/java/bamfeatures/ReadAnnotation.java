package bamfeatures;

import gtf.structs.Exon;
import gtf.structs.Gene;
import gtf.structs.Interval;
import gtf.structs.Transcript;
import gtf.treecollections.IntervalTreeForestManager;
import gtf.treecollections.PCRIndexManager;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import readsimulator.Pair;

import java.util.*;

import static runners.ReadSimulatorRunner.cut;

public class ReadAnnotation {
    public static final char PIPE = '|';
    public static final String MERGED = "MERGED";
    public static final String INTRONIC = "INTRON";
    private static final String SPLIT_INCONSISTENT = "\tsplit-inconsistent:true";
    private static final String MM = "\tmm:";
    private static final String CLIPPING = "\tclipping:";
    private static final String NSPLIT = "\tnsplit:";
    private static final String GCOUNT = "\tgcount:";
    private static final String GDIST = "\tgdist:";
    private static final String ANTISENSE_FALSE = "\tantisense:false";
    private static final String ANTISENSE_TRUE = "\tantisense:true";
    private static final String PCR_INDEX = "\tpcrindex: ";
    private static final char TAB = '\t';
    private static final char COMMA = ',';
    private static final char COLON = ':';
    private TreeSet<Interval> firstRead;
    private TreeSet<Interval> secondRead;
    private boolean transcriptomicProcess = false;
    private final String readID;
    private int splitCount;
    private int clippingSum;
    private int mismatchCount;
    private int geneCount;
    private int pcrIndex;
    private List<Gene> genesThatInclude;
    private TreeSet<Interval> combinedRead;
    private boolean isConsistent;
    // Tree possible result scenarios, but the merged and intronic can be combined
    private List<Pair<Gene, List<Transcript>>> transcriptomicMatches;
    private List<Gene> mergedTranscriptomicMatches;
    private boolean isReadStrandNegative;
    private int alignmentStart;
    private int alignmentEnd;
    private int firstAlignmentStart;
    private int firstAlignmentEnd;
    private int secondAlignmentStart;
    private int secondAlignmentEnd;
    private boolean antisense;
    private int gdist;
    public ReadAnnotation(String readName) {
        this.readID = readName;
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

    public String getReadID() {
        return readID;
    }

    public List<Gene> getGenesThatInclude() {
        return genesThatInclude;
    }

    public void setGenesThatInclude(List<Gene> genesThatInclude) {
        this.genesThatInclude = genesThatInclude;
    }

    public List<Pair<Gene, List<Transcript>>> getTranscriptomicMatches() {
        return transcriptomicMatches;
    }

    public void extractReadAlignmentStartEnd(SAMRecord first, SAMRecord second) {

        // Reads can be flipped in the alignment
        isReadStrandNegative = first.getReadNegativeStrandFlag();
        firstAlignmentStart = first.getAlignmentStart();
        firstAlignmentEnd = first.getAlignmentEnd();
        secondAlignmentStart = second.getAlignmentStart();
        secondAlignmentEnd = second.getAlignmentEnd();

        alignmentStart = Math.min(firstAlignmentStart, secondAlignmentStart);
        alignmentEnd = Math.max(firstAlignmentEnd, secondAlignmentEnd);


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

        isConsistent = checkConsistency(firstRead, secondRead) && checkConsistency(secondRead, firstRead);
        return isConsistent;
    }

    private boolean checkConsistency(TreeSet<Interval> readA, TreeSet<Interval> readB) {
        // Get the splits for readA
        // TODO: don't call up getsplits since when need to calculate them anyway
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

    public void calculateGeneDistance(IntervalTreeForestManager manager) {
        gdist = manager.getDistanceToNearestNeighborGene(this);
    }

    public int getGdist() {
        return gdist;
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

    public void calculateMismatches(SAMRecord first, SAMRecord second) {
// Attribute order NM, nM, then XM
        Integer firstCount = (Integer) first.getAttribute("NM");
        firstCount = (firstCount != null) ? firstCount : (Integer) first.getAttribute("nM");
        firstCount = (firstCount != null) ? firstCount : (Integer) first.getAttribute("XM");
        Integer secondCount = (Integer) second.getAttribute("NM");
        secondCount = (secondCount != null) ? secondCount : (Integer) second.getAttribute("nM");
        secondCount = (secondCount != null) ? secondCount : (Integer) second.getAttribute("XM");
        mismatchCount = firstCount + secondCount;
    }

    public void calculateClipping(SAMRecord first, SAMRecord second) {
        clippingSum = firstAlignmentStart - first.getUnclippedStart();
        clippingSum += first.getUnclippedEnd() - firstAlignmentEnd;
        clippingSum += secondAlignmentStart - second.getUnclippedStart();
        clippingSum += second.getUnclippedEnd() - secondAlignmentEnd;
    }

    public TreeSet<Interval> getCombinedRead() {
        return combinedRead;
    }

    // TODO: merge caluclation of split with inconsistency calculation
    public void calculateSplit() {
        // Split count: the size of the unique implied intron set
        //(fw and rw may imply the same intron(s))
        Interval prev = null;
        splitCount = 0;
        Set<Interval> splits = new HashSet<>();

        // TreeSet to hold combined intervals
        TreeSet<Interval> combinedRead = new TreeSet<>();

        // Process firstRead
        for (Interval interval : firstRead) {
            if (prev != null) {
                splits.add(new Interval(prev.getEnd() + 1, interval.getStart() - 1));
            }
            mergeInterval2(combinedRead, interval);
            prev = interval;
        }

        prev = null;
        // Process secondRead
        for (Interval interval : secondRead) {
            if (prev != null) {
                splits.add(new Interval(prev.getEnd() + 1, interval.getStart() - 1));
            }
            mergeInterval2(combinedRead, interval);
            prev = interval;
        }

        splitCount = splits.size();
        this.combinedRead = combinedRead; // Assuming combinedRead is a member variable
    }

    /**
     * Merges an interval into a TreeSet of intervals, ensuring no overlaps or adjacency.
     */
    private void mergeInterval(TreeSet<Interval> combinedRead, Interval newInterval) {
        Interval lower = combinedRead.floor(newInterval);
        Interval higher = combinedRead.ceiling(newInterval);

        if (lower != null && lower.getEnd() >= newInterval.getStart() - 1) {
            // Merge with the lower interval if overlapping or adjacent
            newInterval = new Interval(lower.getStart(), Math.max(lower.getEnd(), newInterval.getEnd()));
            combinedRead.remove(lower);
        }
        if (higher != null && higher.getStart() <= newInterval.getEnd() + 1) {
            // Merge with the higher interval if overlapping or adjacent
            newInterval = new Interval(newInterval.getStart(), Math.max(higher.getEnd(), newInterval.getEnd()));
            combinedRead.remove(higher);
        }
        // Add the merged interval
        combinedRead.add(newInterval);
    }

    // TODO Check if can combined
    private void mergeInterval2(TreeSet<Interval> combinedRead, Interval newInterval) {
        boolean merged;
        do {
            merged = false;
            Interval lower = combinedRead.floor(newInterval);
            Interval higher = combinedRead.ceiling(newInterval);

            if (lower != null && lower.getEnd() >= newInterval.getStart() - 1) {
                // Merge with the lower interval
                newInterval = new Interval(lower.getStart(), Math.max(lower.getEnd(), newInterval.getEnd()));
                combinedRead.remove(lower);
                merged = true;
            }
            if (higher != null && higher.getStart() <= newInterval.getEnd() + 1) {
                // Merge with the higher interval
                newInterval = new Interval(newInterval.getStart(), Math.max(higher.getEnd(), newInterval.getEnd()));
                combinedRead.remove(higher);
                merged = true;
            }
        } while (merged);
        combinedRead.add(newInterval);
    }

    public void calculatePCRIndex(PCRIndexManager pcrIndex) {
        this.pcrIndex = pcrIndex.getPCRIndex(combinedRead, isReadStrandNegative);
    }

    public int getGeneCount() {
        return geneCount;
    }

    public void processAntisense(IntervalTreeForestManager manager) {
        antisense = manager.isAntisenseBetter(this);
    }

    public boolean isTranscriptomicProcess() {
        return transcriptomicProcess;
    }

    public boolean findTranscriptomicMatches() {
        transcriptomicMatches = new ArrayList<>();
        transcriptomicProcess = true;

        boolean isFound;
        List<Transcript> matchingTranscripts = null;
        for (Gene gene : genesThatInclude) {
            isFound = false;
            // Iterate over each transcript of the gene
            for (Transcript transcript : gene.getTranscripts().values()) {
                TreeSet<Exon> exons = transcript.getExons();
                // TODO: without cutting as this requires a new object
                if (cut(exons, new Interval(firstRead.getFirst().getStart(), firstRead.getLast().getEnd())).equals(firstRead) && cut(exons, new Interval(secondRead.getFirst().getStart(), secondRead.getLast().getEnd())).equals(secondRead)) {
                    if (isFound) {
                        matchingTranscripts.add(transcript);
                    } else {
                        matchingTranscripts = new ArrayList<>();
                        matchingTranscripts.add(transcript);
                        isFound = true;
                    }
                }
            }

            if (isFound) {
                transcriptomicMatches.add(new Pair<>(gene, matchingTranscripts));
            }
        }
        geneCount = transcriptomicMatches.size();
        return !transcriptomicMatches.isEmpty();
    }


    public boolean findMergedTranscriptomicMatches() {
        mergedTranscriptomicMatches = new ArrayList<>();
        // Iterate over all genes
        for (Gene gene : genesThatInclude) {
            // Merge all exons of all transcripts into a single TreeSet<Interval>
            TreeSet<Interval> mergedGeneTranscript = new TreeSet<>();
            for (Transcript transcript : gene.getTranscripts().values()) {
                for (Exon exon : transcript.getExons()) {
                    mergeInterval2(mergedGeneTranscript, exon.getInterval());
                }
            }

            // Check if the combinedRead is fully contained within the mergedGeneTranscript
            if (isReadContainedInMergedTranscript(combinedRead, mergedGeneTranscript)) {
                mergedTranscriptomicMatches.add(gene);
            }
        }

        // Return true if any gene matches
        geneCount = mergedTranscriptomicMatches.size();
        if (geneCount == 0) {
            geneCount = genesThatInclude.size();
        }
        return !mergedTranscriptomicMatches.isEmpty();
    }

    private boolean isReadContainedInMergedTranscript(TreeSet<Interval> readIntervals, TreeSet<Interval> mergedGeneTranscript) {
        Iterator<Interval> readIter = readIntervals.iterator();
        Iterator<Interval> mergedIter = mergedGeneTranscript.iterator();

        Interval currentRead = readIter.hasNext() ? readIter.next() : null;
        Interval currentMerged = mergedIter.hasNext() ? mergedIter.next() : null;

        while (currentRead != null) {
            // Ensure the current merged interval can fully contain the current read interval
            while (currentMerged != null && currentMerged.getEnd() < currentRead.getStart()) {
                currentMerged = mergedIter.hasNext() ? mergedIter.next() : null;
            }

            if (currentMerged == null || currentRead.getStart() < currentMerged.getStart() || currentRead.getEnd() > currentMerged.getEnd()) {
                // No interval in mergedGeneTranscript can fully contain the current read interval
                return false;
            }

            // Move to the next read interval
            currentRead = readIter.hasNext() ? readIter.next() : null;
        }

        // All read intervals were contained
        return true;
    }

    public void findIntronicGenes() {
        geneCount = genesThatInclude.size();
    }

    public int getPcrIndex() {
        return pcrIndex;
    }

    public int getSplitCount() {
        return splitCount;
    }

    public int getClippingSum() {
        return clippingSum;
    }

    public int getMismatchCount() {
        return mismatchCount;
    }

    public boolean isReadStrandNegative() {
        return isReadStrandNegative;
    }

    public int getAlignmentStart() {
        return alignmentStart;
    }

    public int getAlignmentEnd() {
        return alignmentEnd;
    }

    public int getFirstAlignmentStart() {
        return firstAlignmentStart;
    }

    public int getFirstAlignmentEnd() {
        return firstAlignmentEnd;
    }

    public int getSecondAlignmentStart() {
        return secondAlignmentStart;
    }

    public int getSecondAlignmentEnd() {
        return secondAlignmentEnd;
    }

    public List<Gene> getMergedTranscriptomicMatches() {
        return mergedTranscriptomicMatches;
    }

    public String output() {
        StringBuilder sb = new StringBuilder();
        sb.append(readID);
        if (!isConsistent) {
            return sb.append(SPLIT_INCONSISTENT).toString();
        } else {
            sb.append(MM).append(mismatchCount)
                    .append(CLIPPING).append(clippingSum)
                    .append(NSPLIT).append(splitCount)
                    .append(GCOUNT).append(geneCount);
            if (genesThatInclude != null) {
                // transcriptomic
                if (!transcriptomicMatches.isEmpty()) {
                    Iterator<Pair<Gene, List<Transcript>>> it = transcriptomicMatches.iterator();
                    sb.append(TAB);
// First match is handled separately to avoid leading delimiter
                    if (it.hasNext()) {
                        Pair<Gene, List<Transcript>> match = it.next();
                        appendMatch(sb, match);

                        // Process remaining matches
                        while (it.hasNext()) {
                            match = it.next();
                            sb.append(PIPE);
                            appendMatch(sb, match);
                        }
                    }

                } else if (!mergedTranscriptomicMatches.isEmpty()) {
                    // merged YAL005C,protein_coding:MERGED|YAL004W,protein_coding:MERGED
                    outputNonTranscriptomic(sb, mergedTranscriptomicMatches, MERGED);
                } else {
                    // intronic ENSG00000269831,protein_coding:INTRON|ENSG00000230092,pseudogene:INTRON
                    outputNonTranscriptomic(sb, genesThatInclude, INTRONIC);

                }

            } else {
                sb.append(GDIST).append(gdist).append(antisense ? ANTISENSE_TRUE : ANTISENSE_FALSE);
            }
            sb.append(PCR_INDEX).append(pcrIndex);
            return sb.toString();
        }

    }

    private void appendMatch(StringBuilder sb, Pair<Gene, List<Transcript>> match) {
        sb
                .append(match.getFirst().getGeneID())
                .append(COMMA)
                .append(match.getFirst().getSource())
                .append(COLON);

        Iterator<Transcript> transcriptIt = match.getSecond().iterator();
        if (transcriptIt.hasNext()) {
            sb.append(transcriptIt.next().getTranscriptID());
            while (transcriptIt.hasNext()) {
                sb.append(COMMA).append(transcriptIt.next().getTranscriptID());
            }
        }
    }

    private void outputNonTranscriptomic(StringBuilder sb, List<Gene> genesThatInclude, String intronic) {
        Iterator<Gene> it = genesThatInclude.iterator();
        Gene gene = it.next();
        sb.append(TAB).append(gene.getGeneID()).append(COMMA).append(gene.getSource()).append(COLON).append(intronic);
        while (it.hasNext()) {
            gene = it.next();
            sb.append(PIPE).append(gene.getGeneID()).append(COMMA).append(gene.getSource()).append(COLON).append(intronic);
        }
    }
}
