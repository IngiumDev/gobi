package gtf;

import java.util.HashMap;
import java.util.Map;

public class Gene extends AnnotationEntry {
    private String id;
    private Map<String, Transcript> transcripts;
    private Map<String, String> protein_to_transcript;

    public Gene(String seqname, String source, String feature, Interval interval, double score, StrandDirection strand, FrameStarts frame, Map<String, String> attributes) {
        super(seqname, source, feature, interval, score, strand, frame, attributes);
        transcripts = new HashMap<>();
    }

    public Gene(String id) {
        this.id = id;
        transcripts = new HashMap<>();
    }

    public Gene(String id, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        super(seqname, source, strand, attributes);
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

    public void setTranscripts(Map<String, Transcript> transcripts) {
        this.transcripts = transcripts;
    }

    public Transcript getTranscript(String id) {
        return transcripts.get(id);
    }

    public void addTranscript(String id, Transcript transcript) {
        transcripts.put(id, transcript);

    }

    public void processProteins() {
        protein_to_transcript = new HashMap<>();
        for (Transcript transcript : transcripts.values()) {
            if (!transcript.getCds().isEmpty()) {
                String protein_id = transcript.getCds().getFirst().getAttribute("protein_id");
                if (protein_id != null) {
                    protein_to_transcript.put(protein_id, transcript.getId());
                }
            }
        }
    }

    public String convertProteinToTranscript(String protein_id) {
        return protein_to_transcript.get(protein_id);
    }


}
