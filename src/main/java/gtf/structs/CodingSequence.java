package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.Map;

public class CodingSequence extends AnnotationEntry implements Comparable<CodingSequence> {
    private String proteinID;
    private String ccdsID; // Not always

    public CodingSequence(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Attribute attribute) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.proteinID = attribute.getProteinID();
        this.ccdsID = attribute.getCcdsID();
    }

    public String getProteinID() {
        return proteinID;
    }

    public String getCcdsID() {
        return ccdsID;
    }

    @Override
    public int compareTo(CodingSequence o) {
        return Integer.compare(this.getInterval().getStart(), o.getInterval().getStart());
    }
}
