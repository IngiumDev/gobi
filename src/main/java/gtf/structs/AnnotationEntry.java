package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

public abstract class AnnotationEntry {
    // GTF line columns to escalate: seqname, source, strand

    // name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    private String seqname;

    public String getSource() {
        return source;
    }

    public String getFeature() {
        return feature;
    }

    public double getScore() {
        return score;
    }

    public FrameStarts getFrame() {
        return frame;
    }

    // name of the program that generated this feature, or the data source (database or project name)
    private String source;
    // feature type name, e.g. gtf.structs.Gene, Variation, Similarity
    private String feature;
    // gtf.structs.Interval
    private Interval interval;
    // score - A floating point value.
    private double score;
    // strand - defined as + (forward) or - (reverse).
    private StrandDirection strand;
    // frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    private FrameStarts frame;


    public AnnotationEntry(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame) {
        this.seqname = seqname;
        this.source = source;
        this.feature = feature;
        this.interval = interval;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
    }

    public AnnotationEntry(String seqname, String source, StrandDirection strand) {
        this.seqname = seqname;
        this.source = source;
        this.strand = strand;
    }

    public AnnotationEntry() {
    }

    public String getSeqname() {
        return seqname;
    }

    public Interval getInterval() {
        return interval;
    }

    public StrandDirection getStrand() {
        return strand;
    }

    public void overwrite(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame) {
        this.seqname = seqname;
        this.source = source;
        this.feature = feature;
        this.interval = interval;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
    }
}
