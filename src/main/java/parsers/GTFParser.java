package parsers;

import gtf.CodingSequence;
import gtf.GTFAnnotation;
import gtf.structs.Exon;
import gtf.structs.Gene;
import gtf.structs.Interval;
import gtf.structs.Transcript;
import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

import static FileUtils.FileReader.parseAttributes;

public class GTFParser {

    public static final String GENE_ID = "gene_id";
    public static final String TRANSCRIPT_ID = "transcript_id";
    public static final String GENE_NAME = "gene_name";
    public static final String GENE_SOURCE = "gene_source";
    public static final String GENE_BIOTYPE = "gene_biotype";
    public static final String TRANSCRIPT_NAME = "transcript_name";
    public static final String TRANSCRIPT_SOURCE = "transcript_source";
    public final static String[] ATTRIBUTES_TO_ESCALATE_TO_TRANSCRIPT = {GENE_ID, GENE_NAME, GENE_SOURCE, GENE_BIOTYPE, TRANSCRIPT_ID, TRANSCRIPT_NAME, TRANSCRIPT_SOURCE};
    public final static String[] ATTRIBUTES_TO_ESCALATE_TO_EXON = ATTRIBUTES_TO_ESCALATE_TO_TRANSCRIPT;
    public final static String[] ATTRIBUTES_TO_ESCALATE_TO_GENE = {GENE_ID, GENE_NAME, GENE_SOURCE, GENE_BIOTYPE};
    public static final int FEATURE_COL = 2;
    public static final int SEQNAME_COL = 0;
    public static final int SOURCE_COL = 1;
    public static final int STRAND_COL = 6;
    public static final int START_COL = 3;
    public static final int END_COL = 4;
    public static final int SCORE_COL = 5;
    public static final int FRAME_COL = 7;
    public static final int ATTRIBUTE_COL = 8;
    public static final String EXON_NUMBER = "exon_number";
    public static final String PROTEIN_ID = "protein_id";

    // GTF line columns to escalate: seqname, source, strand
    public static GTFAnnotation parseGTF(String gtfFile) {
        GTFAnnotation GTFAnnotation = new GTFAnnotation();
        try (BufferedReader br = new BufferedReader(new FileReader(gtfFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.charAt(1) != '#') {
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
        switch (data[FEATURE_COL]) {
            case "gene" -> processGene(data, GTFAnnotation);
            case "transcript" -> processTranscript(data, GTFAnnotation);
            case "exon" -> processExon(data, GTFAnnotation);
            case "CDS" -> processCDS(data, GTFAnnotation);
        }
    }

    private static void processCDS(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[ATTRIBUTE_COL]);
        String geneId = attributes.get(GENE_ID);
        String transcriptId = attributes.get(TRANSCRIPT_ID);

        Gene gene = getOrCreateGene(GTFAnnotation, geneId, data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);
        Transcript transcript = getOrCreateTranscript(gene, transcriptId, data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);
        int exonNumber = Integer.parseInt(attributes.get(EXON_NUMBER));

        Exon exon = getOrCreateExon(exonNumber, transcript, data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);

        createNewCDS(data, attributes, exon, transcript);


    }

    private static void createNewCDS(String[] data, Map<String, String> attributes, Exon exon, Transcript transcript) {
        CodingSequence cds = new CodingSequence(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        exon.setCds(cds);
        transcript.addCds(cds);
    }

    private static Exon getOrCreateExon(int exonNumber, Transcript transcript, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        Exon exon = transcript.getExonByNumber(exonNumber);
        if (exon == null) {
            exon = new Exon(exonNumber, seqname, source, strand, escalateAttributes(attributes, ATTRIBUTES_TO_ESCALATE_TO_EXON));
            transcript.addExon(exon);
        }
        return exon;
    }

    public static Map<String, String> escalateAttributes(Map<String, String> attributes, String[] keys) {
        Map<String, String> escalated_attributes = new HashMap<>();
        for (String key : keys) {
            if (attributes.containsKey(key)) {
                escalated_attributes.put(key, attributes.get(key));
            }
        }
        return escalated_attributes;
    }

    private static void processExon(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[ATTRIBUTE_COL]);
        String geneId = attributes.get(GENE_ID);
        String transcriptId = attributes.get(TRANSCRIPT_ID);

        Gene gene = getOrCreateGene(GTFAnnotation, geneId, data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);
        Transcript transcript = getOrCreateTranscript(gene, transcriptId, data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);

        int exonNumber = Integer.parseInt(attributes.get(EXON_NUMBER));
        Exon exon = transcript.getExonByNumber(exonNumber);
        if (exon == null) {
            createNewExon(data, attributes, transcript);
        } else {
            overwriteExon(data, exon, attributes);
        }

    }

    private static void overwriteExon(String[] data, Exon exon, Map<String, String> attributes) {
        exon.overwrite(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
    }

    private static void createNewExon(String[] data, Map<String, String> attributes, Transcript transcript) {
        Exon exon;
        exon = new Exon(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        exon.setExonNumber(Integer.parseInt(attributes.get(EXON_NUMBER)));
        transcript.addExon(exon);
    }

    private static Transcript getOrCreateTranscript(Gene gene, String transcriptId, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        Transcript transcript = gene.getTranscript(transcriptId);
        if (transcript == null) {
            transcript = new Transcript(transcriptId, seqname, source, strand, escalateAttributes(attributes, ATTRIBUTES_TO_ESCALATE_TO_TRANSCRIPT));
            gene.addTranscript(transcriptId, transcript);
        }
        return transcript;
    }

    private static void processTranscript(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[ATTRIBUTE_COL]);
        String geneId = attributes.get(GENE_ID);
        String transcriptId = attributes.get(TRANSCRIPT_ID);

        Gene gene = getOrCreateGene(GTFAnnotation, geneId, data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);

        Transcript transcript = gene.getTranscript(transcriptId);
        if (transcript == null) {
            createNewTranscript(data, attributes, transcriptId, gene);
        } else {
            overwriteTranscript(data, transcript, attributes);
        }

    }

    private static void overwriteTranscript(String[] data, Transcript transcript, Map<String, String> attributes) {
        transcript.overwrite(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
    }

    private static void createNewTranscript(String[] data, Map<String, String> attributes, String transcriptId, Gene gene) {
        Transcript transcript;
        transcript = new Transcript(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        transcript.setId(transcriptId);
        gene.addTranscript(transcriptId, transcript);
    }

    private static Gene getOrCreateGene(GTFAnnotation GTFAnnotation, String geneId, String seqname, String source, StrandDirection strand, Map<String, String> attributes) {
        Gene gene = GTFAnnotation.getGenes().get(geneId);
        if (gene == null) {
            gene = new Gene(geneId, seqname, source, strand, escalateAttributes(attributes, ATTRIBUTES_TO_ESCALATE_TO_GENE));
            GTFAnnotation.addGene(gene);
        }
        return gene;
    }

    private static void processGene(String[] data, GTFAnnotation GTFAnnotation) {
        Map<String, String> attributes = parseAttributes(data[ATTRIBUTE_COL]);
        Gene thisGene = GTFAnnotation.getGenes().get(attributes.get(GENE_ID));
        // Check if gene already exists
        if (thisGene == null) {
            createNewGene(data, GTFAnnotation, attributes);
        } else {
            // Overwrite the gtf.structs.AnnotationEntry fields with the new one
            overwriteGene(data, thisGene, attributes);
        }

    }

    private static void overwriteGene(String[] data, Gene thisGene, Map<String, String> attributes) {
        thisGene.overwrite(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
    }

    private static void createNewGene(String[] data, GTFAnnotation GTFAnnotation, Map<String, String> attributes) {
        Gene gene = new Gene(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        gene.setId(gene.getAttribute(GENE_ID));
        GTFAnnotation.addGene(gene);
    }

    public static double parseScore(String data) {
        if (Objects.equals(data, ".")) {
            return SEQNAME_COL;
        } else {
            return Double.parseDouble(data);
        }
    }

    public static StrandDirection parseStrand(String data) {
        return switch (data) {
            case "+" -> StrandDirection.FORWARD;
            case "-" -> StrandDirection.REVERSE;
            case null, default -> throw new IllegalArgumentException("Invalid strand direction: " + data);
        };
    }

    public static FrameStarts parseFrame(String data) {
        return switch (data) {
            case "." -> FrameStarts.NONE;
            case "0" -> FrameStarts.ZERO;
            case "1" -> FrameStarts.ONE;
            case "2" -> FrameStarts.TWO;
            default -> throw new IllegalArgumentException("Invalid frame start: " + data);
        };
    }

}
