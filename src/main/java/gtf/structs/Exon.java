package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

public class Exon extends AnnotationEntry implements Comparable<Exon> {
    private final String exonID; // rarely used
    private final String exonNumber;
    private CodingSequence cds;

    // If we read a EXON line
    public Exon(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Attribute attribute) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.exonNumber = attribute.getExonNumber();
        this.exonID = attribute.getExonID();
    }

    // If we read a CDS line
    public Exon(Attribute attribute) {
        this.exonNumber = attribute.getExonNumber();
        this.exonID = attribute.getExonID();
    }

    public String getExonID() {
        return exonID;
    }

    public String getExonNumber() {
        return exonNumber;
    }

    public CodingSequence getCds() {
        return cds;
    }

    public void setCds(CodingSequence cds) {
        this.cds = cds;
    }

    @Override
    public int compareTo(Exon other) {
        return Integer.compare(this.getInterval().getStart(), other.getInterval().getStart());
    }
}
