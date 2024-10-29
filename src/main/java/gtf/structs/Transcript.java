package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.TreeSet;

public class Transcript extends AnnotationEntry {
    private final String transcriptID;
    private final String transcriptName;
    // Sorted by start position
    private final TreeSet<Exon> exons;
    private final TreeSet<CodingSequence> cds;

    // If we read a transcript line
    public Transcript(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, GTFAttributes GTFAttributes) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.transcriptID = GTFAttributes.getTranscriptID();
        this.transcriptName = GTFAttributes.getTranscriptName();
        this.exons = new TreeSet<>();
        this.cds = new TreeSet<>();
    }

    // If we read an exon or CDS line
    public Transcript(String seqname, String source, StrandDirection strand, GTFAttributes GTFAttributes) {
        super(seqname, source, strand);
        this.transcriptID = GTFAttributes.getTranscriptID();
        this.transcriptName = GTFAttributes.getTranscriptName();
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
