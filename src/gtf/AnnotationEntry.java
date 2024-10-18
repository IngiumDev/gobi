package gtf;

import java.util.Map;

public abstract class AnnotationEntry {
    // name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    private String seqname;
    // name of the program that generated this feature, or the data source (database or project name)
    private String source;
    // feature type name, e.g. gtf.Gene, Variation, Similarity
    private String feature;
    // gtf.Interval
    private Interval interval;
    // score - A floating point value.
    private double score;
    // strand - defined as + (forward) or - (reverse).
    private StrandDirection strand;
    // frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    private FrameStarts frame;
    private Map<String, String> attributes;

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
