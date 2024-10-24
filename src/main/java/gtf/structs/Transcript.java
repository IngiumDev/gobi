package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class Transcript extends AnnotationEntry {
    // TODO: think about where to store cds once
    private String transcriptID;
    private String transcriptName;
    // Sorted by start position
    private TreeSet<Exon> exons;
    private TreeSet<CodingSequence> cds;

    // If we read a transcript line
    public Transcript(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Attribute attribute) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.transcriptID = attribute.getTranscriptID();
        this.transcriptName = attribute.getTranscriptName();
        this.exons = new TreeSet<>();
        this.cds = new TreeSet<>();
    }
    // If we read an exon or CDS line
    public Transcript(Attribute attribute) {
        this.transcriptID = attribute.getTranscriptID();
        this.transcriptName = attribute.getTranscriptName();
        this.exons = new TreeSet<>();
        this.cds = new TreeSet<>();
    }

    public boolean addExon(Exon exon) {
        return exons.add(exon);
    }
    public boolean addCds(CodingSequence cds) {
        return this.cds.add(cds);
    }

    public String getTranscriptName() {
        return transcriptName;
    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public TreeSet<Exon> getExons() {
        return exons;
    }

    public TreeSet<CodingSequence> getCds() {
        return cds;
    }
}
