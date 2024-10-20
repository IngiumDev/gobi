package parsers;

import gtf.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Map;
import java.util.Objects;

import static FileUtils.FileReader.parseAttributes;

public class GTFParser {

    public static Annotation parseGTF(String gtfFile) {
        Annotation annotation = new Annotation();
        try (BufferedReader br = new BufferedReader(new FileReader(gtfFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.charAt(0) != '#') {
                    String[] data = line.split("\t");
                    switch (data[2]) {
                        case "gene" -> processGene(data, annotation);
                        case "transcript" -> processTranscript(data, annotation);
                        case "exon" -> processExon(data, annotation);
                        case "CDS" -> processCDS(data, annotation);
                    }

                }

            }
            // Process introns
            annotation.getGenes().values().parallelStream().forEach(gene -> {
                gene.getTranscripts().values().parallelStream().forEach(Transcript::processIntrons);
            });

        } catch (Exception e) {
            e.printStackTrace();
        }
        return annotation;
    }

    private static void processCDS(String[] data, Annotation annotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        String geneId = attributes.get("gene_id");
        String transcriptId = attributes.get("transcript_id");

        Gene gene = getOrCreateGene(annotation, geneId);
        Transcript transcript = getOrCreateTranscript(gene, transcriptId);
        int exonNumber = Integer.parseInt(attributes.get("exon_number"));
        Exon exon = transcript.getExonByNumber(exonNumber);
        exon = getOrCreateExon(exon, exonNumber, transcript);

        CodingSequence cds = new CodingSequence(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        exon.setCds(cds);
        transcript.addCds(cds);


    }

    private static Exon getOrCreateExon(Exon exon, int exonNumber, Transcript transcript) {
        if (exon == null) {
            exon = new Exon();
            exon.setExonNumber(exonNumber);
            transcript.addExon(exon);
        }
        return exon;
    }

    private static void processExon(String[] data, Annotation annotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        String geneId = attributes.get("gene_id");
        String transcriptId = attributes.get("transcript_id");

        Gene gene = getOrCreateGene(annotation, geneId);
        Transcript transcript = getOrCreateTranscript(gene, transcriptId);

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

    private static Transcript getOrCreateTranscript(Gene gene, String transcriptId) {
        Transcript transcript = gene.getTranscript(transcriptId);
        if (transcript == null) {
            transcript = new Transcript(transcriptId);
            gene.addTranscript(transcriptId, transcript);
        }
        return transcript;
    }

    private static void processTranscript(String[] data, Annotation annotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        String geneId = attributes.get("gene_id");
        String transcriptId = attributes.get("transcript_id");

        Gene gene = getOrCreateGene(annotation, geneId);

        Transcript transcript = gene.getTranscript(transcriptId);
        if (transcript == null) {
            transcript = new Transcript(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
            transcript.setId(transcriptId);
            gene.addTranscript(transcriptId, transcript);
        } else {
            transcript.overwrite(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        }

    }

    private static Gene getOrCreateGene(Annotation annotation, String geneId) {
        Gene gene = annotation.getGenes().get(geneId);
        if (gene == null) {
            gene = new Gene(geneId);
            annotation.addGene(gene);
        }
        return gene;
    }

    private static void processGene(String[] data, Annotation annotation) {
        Map<String, String> attributes = parseAttributes(data[8]);
        Gene thisGene = annotation.getGenes().get(attributes.get("gene_id"));
        // Check if gene already exists
        if (thisGene != null) {
            // Overwrite the gtf.AnnotationEntry fields with the new one
            thisGene.overwrite(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
        } else {
            Gene gene = new Gene(data[0], data[1], data[2], new Interval(Integer.parseInt(data[3]), Integer.parseInt(data[4])), getScore(data[5]), getStrand(data[6]), getFrame(data[7]), attributes);
            gene.setId(gene.getAttribute("gene_id"));
            annotation.addGene(gene);
        }

    }

    public static double getScore(String data) {
        if (Objects.equals(data, ".")) {
            return 0;
        } else {
            return Double.valueOf(data);
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
