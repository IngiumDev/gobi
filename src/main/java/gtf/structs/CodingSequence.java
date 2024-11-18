package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

public class CodingSequence extends AnnotationEntry implements Comparable<CodingSequence> {
    private final String proteinID;
    private final String ccdsID; // Not always
    private final String exonNumber;

    public CodingSequence(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, GTFAttributes GTFAttributes) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.proteinID = GTFAttributes.getProteinID();
        this.ccdsID = GTFAttributes.getCcdsID();
        this.exonNumber = GTFAttributes.getExonNumber();
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
