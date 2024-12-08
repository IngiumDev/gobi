package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

public class Exon extends AnnotationEntry implements Comparable<Exon> {
    private final String exonID; // rarely used
    private final String exonNumber;

    // If we read a EXON line
    public Exon(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, GTFAttributes GTFAttributes) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.exonNumber = GTFAttributes.getExonNumber();
        this.exonID = GTFAttributes.getExonID();
    }

    @Override
    public int compareTo(Exon other) {
        int startComparison = Integer.compare(this.getInterval().getStart(), other.getInterval().getStart());
        if (startComparison != 0) {
            return startComparison;
        }
        return Integer.compare(this.getInterval().getEnd(), other.getInterval().getEnd());
    }

    @Override
    public String toString() {
        return "Exon " + exonNumber + " " + getInterval().getStart() + "-" + getInterval().getEnd() + " Length: " + getInterval().getLength();
    }
}
