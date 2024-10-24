package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.*;

public class Gene extends AnnotationEntry {
    private String geneID;
    private String geneName;
    private Map<String, Transcript> transcripts;
    // Only when process introns is run
    private Set<Interval> introns;

    //If we read a gene line
    public Gene(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Attribute attribute) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.geneID = attribute.getGeneID();
        this.geneName = attribute.getGeneName();
        this.transcripts = new HashMap<>();
    }
    //If we read a transcript/exon/cds line
    public Gene(Attribute attribute) {
        this.geneID = attribute.getGeneID();
        this.geneName = attribute.getGeneName();
        this.transcripts = new HashMap<>();
    }
    public void addTranscript(Transcript transcript) {
        transcripts.put(transcript.getTranscriptID(), transcript);
    }

    public void processIntrons() {
        introns = new TreeSet<>();
        for (Transcript transcript : transcripts.values()) {
            CodingSequence previousCds = null;
            for (CodingSequence cds : transcript.getCds()) {
                if (previousCds != null) {
                    int intronStart = previousCds.getInterval().getEnd() + 1;
                    int intronEnd = cds.getInterval().getStart() - 1;
                    if (intronStart <= intronEnd) {
                        introns.add(new Interval(intronStart, intronEnd));
                    }
                }
                previousCds = cds;
            }
        }
    }
}
