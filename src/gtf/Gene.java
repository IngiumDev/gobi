package gtf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Gene extends AnnotationEntry {
    private String id;
    private Map<String, Transcript> transcripts;
    //private ArrayList<Proteins>

    public Gene(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        super(seqname, source, feature, interval, score, strand, frame, attributes);
        transcripts = new HashMap<>();
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public Map<String, Transcript> getTranscripts() {
        return transcripts;
    }
    public Transcript getTranscript(String id) {
        return transcripts.get(id);
    }

    public void setTranscripts(Map<String, Transcript> transcripts) {
        this.transcripts = transcripts;
    }

    public Gene(String id) {
        this.id = id;
        transcripts = new HashMap<>();
    }

    public void addTranscript(String id, Transcript transcript) {
        transcripts.put(id, transcript);

    }
}
