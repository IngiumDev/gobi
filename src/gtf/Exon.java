package gtf;

import java.util.Map;

public class Exon extends AnnotationEntry implements Comparable<Exon> {
    private int exonNumber;
    private CodingSequence cds;

    public Exon() {

    }

    // TODO: exons only have one CDS, so we can remove the TreeSet and use a single CodingSequence object, move list of CDS to Transcript
    @Override
    public int compareTo(Exon other) {
        return Integer.compare(this.getInterval().getStart(), other.getInterval().getStart());
    }



    public int getExonNumber() {
        return exonNumber;
    }

    public void setExonNumber(int exonNumber) {
        this.exonNumber = exonNumber;
    }

    public Exon(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        super(seqname, source, feature, interval, score, strand, frame, attributes);

    }

    public CodingSequence getCds() {
        return cds;
    }

    public void setCds(CodingSequence cds) {
        this.cds = cds;
    }
}
