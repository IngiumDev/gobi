import java.util.HashMap;

public abstract class AnnotationEntry {
    // name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    private String Seqname;
    // name of the program that generated this feature, or the data source (database or project name)
    private String source;
    // feature type name, e.g. Gene, Variation, Similarity
    private String feature;
    // Interval
    private Interval interval;
    // score - A floating point value.
    private double score;
    // strand - defined as + (forward) or - (reverse).
    private StrandDirection strand;
    // frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    private FrameStarts frame;
    private HashMap<String, String> attributes = new HashMap<String, String>();

    public AnnotationEntry(HashMap<String, String> attributes, String feature, FrameStarts frame, Interval interval, double score, String seqname, String source, StrandDirection strand) {
        this.attributes = attributes;
        this.feature = feature;
        this.frame = frame;
        this.interval = interval;
        this.score = score;
        Seqname = seqname;
        this.source = source;
        this.strand = strand;
    }

    public String getAttribute(String key) {
        return attributes.get(key);
    }

}
