package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.Map;

public abstract class AnnotationEntry {
    // TODO decouple strand from annotation entry and move it to Gene
    // TODO Migrate away from hashmap for attributes
    // GTF line columns to escalate: seqname, source, strand

    // name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    private String seqname;
    // name of the program that generated this feature, or the data source (database or project name)
    private String source;
    // feature type name, e.g. gtf.structs.Gene, Variation, Similarity
    private String feature;


    public AnnotationEntry(String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        this.seqname = seqname;
        this.source = source;
        this.strand = strand;
        this.attributes = attributes;
    }

    // gtf.structs.Interval
    private Interval interval;
    // score - A floating point value.
    private double score;
    // strand - defined as + (forward) or - (reverse).
    private StrandDirection strand;
    // frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    private FrameStarts frame;
    private Map<String, String> attributes;

    public String getFeature() {
        return feature;
    }

    public void setFeature(String feature) {
        this.feature = feature;
    }

    public FrameStarts getFrame() {
        return frame;
    }

    public void setFrame(FrameStarts frame) {
        this.frame = frame;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public String getSeqname() {
        return seqname;
    }

    public void setSeqname(String seqname) {
        this.seqname = seqname;
    }

    public String getSource() {
        return source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public StrandDirection getStrand() {
        return strand;
    }

    public void setStrand(StrandDirection strand) {
        this.strand = strand;
    }

    public Map<String, String> getAttributes() {
        return attributes;
    }

    public void setAttributes(Map<String, String> attributes) {
        this.attributes = attributes;
    }

    public AnnotationEntry(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        this.seqname = seqname;
        this.source = source;
        this.feature = feature;
        this.interval = interval;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
        this.attributes = attributes;
    }


    public String getAttribute(String key) {
        return attributes.get(key);
    }

    public AnnotationEntry() {
    }

    public Interval getInterval() {
        return interval;
    }

    public void setInterval(Interval interval) {
        this.interval = interval;
    }

    public void overwrite(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        this.seqname = seqname;
        this.source = source;
        this.feature = feature;
        this.interval = interval;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
        this.attributes = attributes;
    }
}
