package parsers;

import gtf.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

import static FileUtils.FileReader.parseAttributes;

public class GTFParser {

    public final static String[] ATTRIBUTES_TO_ESCALATE_TO_TRANSCRIPT = {"gene_id", "gene_name", "gene_source", "gene_biotype", "transcript_id", "transcript_name", "transcript_source"};
    public final static String[] ATTRIBUTES_TO_ESCALATE_TO_GENE = {"gene_id", "gene_name", "gene_source", "gene_biotype"};
    public final static String[] ATTRIBUTES_TO_ESCALATE_TO_EXON = ATTRIBUTES_TO_ESCALATE_TO_TRANSCRIPT;
    // GTF line columns to escalate: seqname, source, strand
    public static GTFAnnotation parseGTF(String gtfFile) {
        GTFAnnotation GTFAnnotation = new GTFAnnotation();
        try (BufferedReader br = new BufferedReader(new FileReader(gtfFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.charAt(0) != '#') {
                    processGTFLine(line, GTFAnnotation);
                }
            }

            // Process introns
            GTFAnnotation.getGenes().values().parallelStream().forEach(gene -> gene.getTranscripts().values().parallelStream().forEach(Transcript::processIntrons));
            GTFAnnotation.getGenes().values().parallelStream().forEach(Gene::processProteins);

        } catch (Exception e) {
            e.printStackTrace();
        }
        return GTFAnnotation;
    }

    private static void processGTFLine(String line, GTFAnnotation GTFAnnotation) {
        String[] data = line.split("\t");
        switch (data[2]) {
            case "gene" -> processGene(data, GTFAnnotation);
            case "transcript" -> processTranscript(data, GTFAnnotation);
            case "exon" -> processExon(data, GTFAnnotation);
            case "CDS" -> processCDS(data, GTFAnnotation);
        }
    }

    private static void processCDS(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        String geneId = attributes.get("gene_id");
        String transcriptId = attributes.get("transcript_id");

        Gene gene = getOrCreateGene(GTFAnnotation, geneId,data[0], data[1], getStrand(data[6]), attributes);
        Transcript transcript = getOrCreateTranscript(gene, transcriptId, data[0], data[1], getStrand(data[6]), attributes);
        int exonNumber = Integer.parseInt(attributes.get("exon_number"));

        Exon exon = getOrCreateExon(exonNumber, transcript, data[0], data[1], getStrand(data[6]), attributes);

        CodingSequence cds = new CodingSequence(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        exon.setCds(cds);
        transcript.addCds(cds);


    }

    private static Exon getOrCreateExon( int exonNumber, Transcript transcript, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        Exon exon = transcript.getExonByNumber(exonNumber);
        if (exon == null) {
            Map<String, String> escalated_attributes = new HashMap<>();
            for (String key : ATTRIBUTES_TO_ESCALATE_TO_EXON) {
                if (attributes.containsKey(key)) {
                    escalated_attributes.put(key, attributes.get(key));
                }
            }
            exon = new Exon(exonNumber, seqname, source, strand, escalated_attributes);
            transcript.addExon(exon);
        }
        return exon;
    }

    private static void processExon(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        String geneId = attributes.get("gene_id");
        String transcriptId = attributes.get("transcript_id");

        Gene gene = getOrCreateGene(GTFAnnotation, geneId,data[0], data[1], getStrand(data[6]), attributes);
        Transcript transcript = getOrCreateTranscript(gene, transcriptId, data[0], data[1], getStrand(data[6]), attributes);

        int exonNumber = Integer.parseInt(attributes.get("exon_number"));
        Exon exon = transcript.getExonByNumber(exonNumber);
        if (exon == null) {
            exon = new Exon(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        } else {
            exon.overwrite(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        }


        // Possible NPE here if exon_number is not present in the attributes
        exon.setExonNumber(Integer.parseInt(attributes.get("exon_number")));
        transcript.addExon(exon);

    }

    private static Transcript getOrCreateTranscript(Gene gene, String transcriptId, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        Transcript transcript = gene.getTranscript(transcriptId);
        if (transcript == null) {
            HashMap<String, String> escalated_attributes = new HashMap<>();
            for (String key : ATTRIBUTES_TO_ESCALATE_TO_TRANSCRIPT) {
                if (attributes.containsKey(key)) {
                    escalated_attributes.put(key, attributes.get(key));
                }
            }
            transcript = new Transcript(transcriptId, seqname, source, strand, escalated_attributes);
            gene.addTranscript(transcriptId, transcript);
        }
        return transcript;
    }

    private static void processTranscript(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        String geneId = attributes.get("gene_id");
        String transcriptId = attributes.get("transcript_id");

        Gene gene = getOrCreateGene(GTFAnnotation, geneId, data[0], data[1],getStrand(data[6]), attributes);

        Transcript transcript = gene.getTranscript(transcriptId);
        if (transcript == null) {
            transcript = new Transcript(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
            transcript.setId(transcriptId);
            gene.addTranscript(transcriptId, transcript);
        } else {
            transcript.overwrite(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        }

    }

    private static Gene getOrCreateGene(GTFAnnotation GTFAnnotation, String geneId, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        Gene gene = GTFAnnotation.getGenes().get(geneId);
        if (gene == null) {
            // Only copy the attributes in ATTRIBUTES_TO_ESCALATE_TO_GENE
            HashMap<String, String> escalated_attributes = new HashMap<>();
            for (String key : ATTRIBUTES_TO_ESCALATE_TO_GENE) {
                if (attributes.containsKey(key)) {
                    escalated_attributes.put(key, attributes.get(key));
                }
            }
            gene = new Gene(geneId, seqname, source, strand, escalated_attributes);
            GTFAnnotation.addGene(gene);
        }
        return gene;
    }

    private static void processGene(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        Gene thisGene = GTFAnnotation.getGenes().get(attributes.get("gene_id"));
        // Check if gene already exists
        if (thisGene != null) {
            // Overwrite the gtf.AnnotationEntry fields with the new one
            thisGene.overwrite(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        } else {
            Gene gene = new Gene(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
            gene.setId(gene.getAttribute("gene_id"));
            GTFAnnotation.addGene(gene);
        }

    }

    public static double getScore(String data) {
        if (Objects.equals(data, ".")) {
            return 0;
        } else {
            return Double.parseDouble(data);
        }
    }

    public static StrandDirection getStrand(String data) {
        return switch (data) {
            case "+" -> StrandDirection.FORWARD;
            case "-" -> StrandDirection.REVERSE;
            case null, default -> throw new IllegalArgumentException("Invalid strand direction: " + data);
        };
    }

    public static FrameStarts getFrame(String data) {
        return switch (data) {
            case "." -> FrameStarts.NONE;
            case "0" -> FrameStarts.ZERO;
            case "1" -> FrameStarts.ONE;
            case "2" -> FrameStarts.TWO;
            default -> throw new IllegalArgumentException("Invalid frame start: " + data);
        };
    }

}
