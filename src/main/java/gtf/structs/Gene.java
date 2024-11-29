package gtf.structs;

import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class Gene extends AnnotationEntry implements augmentedTree.Interval {
    private final String geneID;
    private final String geneName;
    private final Map<String, Transcript> transcripts;
    // Only when process introns is run
    private Set<Interval> introns;

    //If we read a gene line
    public Gene(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, GTFAttributes GTFAttributes) {
        super(seqname, source, feature, interval, score, strand, frame);
        this.geneID = GTFAttributes.getGeneID();
        this.geneName = GTFAttributes.getGeneName();
        this.transcripts = new HashMap<>();

    }

    //If we read a transcript/exon/cds line
    public Gene(String seqname, String source, StrandDirection strand, GTFAttributes GTFAttributes) {
        super(seqname, source, strand);
        this.geneID = GTFAttributes.getGeneID();
        this.geneName = GTFAttributes.getGeneName();
        this.transcripts = new HashMap<>();
        System.out.println("Warning, gene created with possibly wrong biotype");
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

    public Set<Interval> getIntrons() {
        return introns;
    }

    public Map<String, Transcript> getTranscripts() {
        return transcripts;
    }

    public String getGeneName() {
        return geneName;
    }

    public String getGeneID() {
        return geneID;
    }

    public Transcript getTranscript(String transcriptID) {
        return transcripts.get(transcriptID);
    }

    /**
     * Start position (zero-based inclusive)
     *
     * @return
     */
    @Override
    public int getStart() {
        return this.getInterval().getStart();
    }

    /**
     * End position (zero-based inclusive)
     *
     * @return
     */
    @Override
    public int getStop() {
        return this.getInterval().getEnd();
    }


}
