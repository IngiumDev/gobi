package parsers;

import gtf.GTFAnnotation;
import gtf.structs.*;
import gtf.types.FrameStarts;
import gtf.types.StrandDirection;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Objects;

import static gtf.structs.GTFAttributes.parseAttributes;

public class GTFParser {


    public static final int FEATURE_COL = 2;
    public static final int SEQNAME_COL = 0;
    public static final int SOURCE_COL = 1;
    public static final int STRAND_COL = 6;
    public static final int START_COL = 3;
    public static final int END_COL = 4;
    public static final int SCORE_COL = 5;
    public static final int FRAME_COL = 7;
    public static final int ATTRIBUTE_COL = 8;


    // TODO Logging for errors logback slf4j ch. qos
    // TODO: Remove unnecessary comments/methods
    // TODO: Remove intron parsing
    // GTF line columns to escalate: seqname, source, strand
    public static GTFAnnotation parseGTF(String gtfFile) {
        long startTime = System.currentTimeMillis();
        GTFAnnotation GTFAnnotation = new GTFAnnotation();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(gtfFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.charAt(0) != '#') {
                    processGTFLine(line, GTFAnnotation);
                }
            }
            GTFTimer.setGtfParseTime(System.currentTimeMillis() - startTime);
            System.out.println("LOG: Total time to parse GTF: " + GTFTimer.getGtfParseTime() + " ms");
            // Process introns
            startTime = System.currentTimeMillis();
            GTFAnnotation.getGenes().values().parallelStream().forEach(Gene::processIntrons);

            GTFTimer.setIntronProcessTime(System.currentTimeMillis() - startTime);
            System.out.println("LOG: Total time to process introns: " + GTFTimer.getIntronProcessTime() + " ms");
        } catch (IOException e) {
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

    private static void processCDS(String[] data, GTFAnnotation gtfAnnotation) {
        GTFAttributes attributes = parseAttributes(data[ATTRIBUTE_COL]);
        Gene gene = getOrCreateGene(gtfAnnotation, data, attributes);
        Transcript transcript = getOrCreateTranscript(gene, data, attributes);
        // Assumption: CDS comes after exon, cds does not belong to exon
        CodingSequence cds = new CodingSequence(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        transcript.addCds(cds);
    }

    private static void processExon(String[] data, GTFAnnotation gtfAnnotation) {
        GTFAttributes attributes = parseAttributes(data[ATTRIBUTE_COL]);
        Gene gene = getOrCreateGene(gtfAnnotation, data, attributes);
        Transcript transcript = getOrCreateTranscript(gene, data, attributes);
        // Assumption: CDS comes after exon
        Exon exon = new Exon(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        transcript.addExon(exon);
    }

    private static Transcript getOrCreateTranscript(Gene gene, String[] data, GTFAttributes attributes) {
        Transcript transcript = gene.getTranscript(attributes.getTranscriptID());
        if (transcript == null) {
            transcript = new Transcript(data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);
            gene.addTranscript(transcript);
        }
        return transcript;
    }

    private static void processTranscript(String[] data, GTFAnnotation gtfAnnotation) {
        GTFAttributes attributes = parseAttributes(data[ATTRIBUTE_COL]);

        Gene gene = getOrCreateGene(gtfAnnotation, data, attributes);
        Transcript transcript = gene.getTranscript(attributes.getTranscriptID());
        if (transcript == null) {
            createNewTranscript(data, gtfAnnotation, attributes, gene);
        } else {
            overwriteTranscript(data, transcript);
        }
    }

    private static void overwriteTranscript(String[] data, Transcript transcript) {
        transcript.overwrite(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]));
    }

    private static void createNewTranscript(String[] data, GTFAnnotation gtfAnnotation, GTFAttributes attributes, Gene gene) {
        gene.addTranscript(new Transcript(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes));
    }

    private static Gene getOrCreateGene(GTFAnnotation gtfAnnotation, String[] data, GTFAttributes attributes) {
        Gene gene = gtfAnnotation.getGenes().get(attributes.getGeneID());
        if (gene == null) {
            gene = new Gene(data[SEQNAME_COL], data[SOURCE_COL], parseStrand(data[STRAND_COL]), attributes);
            gtfAnnotation.addGene(gene);
        }
        return gene;
    }

    private static void processGene(String[] data, GTFAnnotation gtfAnnotation) {
        GTFAttributes attributes = parseAttributes(data[ATTRIBUTE_COL]);
        Gene thisGene = gtfAnnotation.getGenes().get(attributes.getGeneID());
        // Check if gene already exists
        if (thisGene == null) {
            createNewGene(data, gtfAnnotation, attributes);
        } else {
            // Overwrite the gtf.structs.AnnotationEntry fields with the new one
            overwriteGene(data, thisGene);
        }
    }

    private static void createNewGene(String[] data, GTFAnnotation gtfAnnotation, GTFAttributes attributes) {
        Gene gene = new Gene(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]), attributes);
        gtfAnnotation.addGene(gene);
    }

    private static void overwriteGene(String[] data, Gene thisGene) {
        thisGene.overwrite(data[SEQNAME_COL], data[SOURCE_COL], data[FEATURE_COL], new Interval(Integer.parseInt(data[START_COL]), Integer.parseInt(data[END_COL])), parseScore(data[SCORE_COL]), parseStrand(data[STRAND_COL]), parseFrame(data[FRAME_COL]));
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
