package gtf;

import java.util.Map;
import java.util.TreeSet;

public class Exon extends AnnotationEntry implements Comparable<Exon> {
    private int exonNumber;
    private TreeSet<CodingSequence> cds;

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
        cds = new TreeSet<>();
    }
    public boolean addCDS(CodingSequence cds) {
        return this.cds.add(cds);
    }

    public TreeSet<CodingSequence> getCds() {
        return cds;
    }

    public void setCds(TreeSet<CodingSequence> cds) {
        this.cds = cds;
    }
}
