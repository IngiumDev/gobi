package gtf;

import java.util.Map;
import java.util.TreeSet;

public class Transcript extends AnnotationEntry {
    private String id;
    // Sorted by start position
    private TreeSet<Exon> exons;
    private TreeSet<CodingSequence> cds;

    public Transcript(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        super(seqname, source, feature, interval, score, strand, frame, attributes);
        exons = new TreeSet<>();
        cds = new TreeSet<>();
    }

    public Transcript(String id) {
        this.id = id;
        exons = new TreeSet<>();
        cds = new TreeSet<>();
    }

    public TreeSet<Exon> getExons() {
        return exons;
    }

    public void setExons(TreeSet<Exon> exons) {
        this.exons = exons;
    }

    public boolean addExon(Exon exon) {
        return exons.add(exon);
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public boolean hasExonWithNumber(int exonNumber) {
        for (Exon exon : exons) {
            if (exon.getExonNumber() == exonNumber) {
                return true;
            }
        }
        return false;
    }

    public Exon getExonByNumber(int exonNumber) {
        for (Exon exon : exons) {
            if (exon.getExonNumber() == exonNumber) {
                return exon;
            }
        }
        return null;
    }

    public boolean addCds(CodingSequence cds) {
        return this.cds.add(cds);
    }

}
