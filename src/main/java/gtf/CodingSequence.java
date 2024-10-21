package gtf;

import gtf.structs.AnnotationEntry;
import gtf.structs.Interval;
import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.Map;

public class CodingSequence extends AnnotationEntry implements Comparable<CodingSequence> {
    public CodingSequence(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        super(seqname, source, feature, interval, score, strand, frame, attributes);
    }

    @Override
    public int compareTo(CodingSequence o) {
        return Integer.compare(this.getInterval().getStart(), o.getInterval().getStart());
    }
}
