package tests;

import bamfeatures.ReadAnnotation;
import bamfeatures.ReadAnnotator;
import gtf.structs.Gene;
import gtf.structs.Transcript;
import gtf.types.StrandDirection;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.junit.jupiter.api.Test;
import readsimulator.Pair;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static bamfeatures.ReadAnnotation.INTRONIC;
import static bamfeatures.ReadAnnotation.MERGED;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class ReadAnnotatorTest {

    public static final String LINE_SEPARATOR = "\n----------------------------------";

    private static List<ReferenceEntry> readReferenceTable(File refTable) throws IOException {
        List<ReferenceEntry> referenceEntries = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(refTable))) {
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                // Skip the header row


                String[] parts = line.split("\t");
                if (parts.length < 4) continue; // Skip invalid rows

                String bam = parts[0];
                String gtf = parts[1];
                Boolean strandness = parts[2].isEmpty() ? null : Boolean.parseBoolean(parts[2]);
                String referenceSolution = parts[3];

                referenceEntries.add(new ReferenceEntry(bam, gtf, strandness, referenceSolution));
            }
        }
        return referenceEntries;
    }

    private static void checkIfGeneralStatsMatch(String[] outputParts, String[] parts, String line, String output) {
        assertEquals(outputParts[1], parts[1], "Mismatch does not match the reference solution\n" +
                "Reference solution: " + line + "\n" +
                "Actual solution: " + output + LINE_SEPARATOR);
        assertEquals(outputParts[2], parts[2], "Clipping does not match the reference solution\n" +
                "Reference solution: " + line + "\n" +
                "Actual solution: " + output + LINE_SEPARATOR);
        assertEquals(outputParts[3], parts[3], "Nsplit does not match the reference solution\n" +
                "Reference solution: " + line + "\n" +
                "Actual solution: " + output + LINE_SEPARATOR);
        assertEquals(outputParts[4], parts[4], "Gcount does not match the reference solution\n" +
                "Reference solution: " + line + "\n" +
                "Actual solution: " + output + LINE_SEPARATOR);
        assertEquals(outputParts[6], parts[6], "Pcrindex does not match the reference solution\n" +
                "Reference solution: " + line + "\n" +
                "Actual solution: " + output + LINE_SEPARATOR);
    }

    @Test
    public void testReadAnnotator() throws IOException {
        Path bamPath = Paths.get("data", "bamfeatures");
        Path gtfPath = Paths.get("data", "gtf");
        File refTable = new File(String.valueOf(Paths.get("data", "bamfeatures", "ref.table")));

        List<ReferenceEntry> referenceEntries = readReferenceTable(refTable);

        // Iterate and print the entries
        // Preprocess gtfs so that we don't have to do it for each read
        referenceEntries.parallelStream().forEach(x -> checkCorrectness(x, bamPath, gtfPath));

    }

    private void checkCorrectness(ReferenceEntry x, Path bamPath, Path gtfPath) {
        SamReader samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamPath.toString(), x.bam));
        StrandDirection strandSpecificity = x.strandness == null ? StrandDirection.UNSPECIFIED : x.strandness ? StrandDirection.FORWARD : StrandDirection.REVERSE;
        ReadAnnotator readAnnotator = new ReadAnnotator.Builder()
                .setSamReader(samReader)
                .setGtfFile(new File(gtfPath.toString(), x.gtf))
                .setStrandSpecificity(strandSpecificity)
                .build();
        long startTime = System.currentTimeMillis();
        List<ReadAnnotation> results = readAnnotator.annotateAndReturnReads();
        System.out.println("Annotation time\t" + (System.currentTimeMillis() - startTime));
        // CAUTION: there are null results (if it was skipped)
        // Add the non null to a Map by the read name, if there are multiple results for the same read name, error
        Map<String, ReadAnnotation> readAnnotationMap = new HashMap<>();
        for (ReadAnnotation annotation : results) {
            if (annotation != null) {
                String readName = annotation.getReadID();
                if (readAnnotationMap.containsKey(readName)) {
                    throw new IllegalStateException("Duplicate read name found: " + readName);
                }
                readAnnotationMap.put(readName, annotation);
            }
        }
        // Now we go through the reference solution and check if the results are correct
        try (BufferedReader br = new BufferedReader(new FileReader(new File(bamPath.toString(), x.referenceSolution)))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split("\t");
                String readName = parts[0];
                ReadAnnotation readAnnotation = readAnnotationMap.get(readName);
                // Remove
                readAnnotationMap.remove(readName);
                if (readAnnotation == null) {
                    throw new IllegalStateException("Read name not found in the results: " + readName);
                }
                if (!readAnnotation.isTranscriptomicProcess()) {
                    // Check if the whole line is equal
                    assertEquals(line, readAnnotation.output(), "Read annotation does not match the reference solution\n" +
                            "Reference solution: " + line + "\n" +
                            "Actual solution: " + readAnnotation.output() + LINE_SEPARATOR);

                } else {
                    // Since transcriptomic process, the output order could be different, however, if we only have one it should still match
                    if (!readAnnotation.getTranscriptomicMatches().isEmpty()) {
                        if (readAnnotation.getTranscriptomicMatches().size() == 1 && readAnnotation.getTranscriptomicMatches().getFirst().getSecond().size() == 1) {
                            assertEquals(line, readAnnotation.output(), "Read annotation does not match the reference solution\n" +
                                    "Reference solution: " + line + "\n" +
                                    "Actual solution: " + readAnnotation.output() + LINE_SEPARATOR);
                        } else {
                            // split the output by \t
                            String output = readAnnotation.output();
                            String[] outputParts = output.split("\t");
                            // Check if cols 1=mismatch,2=clipping,3=nsplit,4=gcount,6=pcrindex
                            checkIfGeneralStatsMatch(outputParts, parts, line, output);
                            // Now we need to check if the gene matches
                            String[] referenceGenes = parts[5].split("\\|");
                            String[] actualGenes = outputParts[5].split("\\|");
                            assertEquals(referenceGenes.length, actualGenes.length, "Number of genes does not match the reference solution\n" +
                                    "Reference solution: " + line + "\n" +
                                    "Actual solution: " + output + LINE_SEPARATOR);
                            for (String refGen : referenceGenes) {
                                String[] refGeneParts = refGen.split(":");
                                String[] refGeneStats = refGeneParts[0].split(",");
                                String readGene = refGeneStats[0];
                                String bioType = refGeneStats[1];
                                String[] refTranscripts = refGeneParts[1].split(",");
                                boolean isGeneFound = false;
                                for (Pair<Gene, List<Transcript>> result : readAnnotation.getTranscriptomicMatches()) {
                                    if (result.getFirst().getGeneID().equals(readGene)) {
                                        isGeneFound = true;
                                        assertEquals(bioType, result.getFirst().getSource(), "Biotype does not match the reference solution\n" +
                                                "Reference solution: " + line + "\n" +
                                                "Actual solution: " + output + LINE_SEPARATOR);
                                        // check if # of transcripts match
                                        assertEquals(refTranscripts.length, result.getSecond().size(), "Number of transcripts does not match the reference solution\n" +
                                                "Reference solution: " + line + "\n" +
                                                "Actual solution: " + output + LINE_SEPARATOR);
                                        for (Transcript t : result.getSecond()) {
                                            boolean isTranscriptFound = false;
                                            for (String refTranscript : refTranscripts) {
                                                if (t.getTranscriptID().equals(refTranscript)) {
                                                    isTranscriptFound = true;
                                                    break;
                                                }
                                            }
                                            if (!isTranscriptFound) {
                                                throw new IllegalStateException("Transcript not found in the results:\n" +
                                                        "Reference solution: " + line + "\n" +
                                                        "Actual solution: " + output + LINE_SEPARATOR);
                                            }
                                        }
                                    }
                                }
                                if (!isGeneFound) {
                                    throw new IllegalStateException("Gene not found in the results:\n" +
                                            "Reference solution: " + line + "\n" +
                                            "Actual solution: " + output + LINE_SEPARATOR);
                                }
                            }

                        }
                    } else if (!readAnnotation.getMergedTranscriptomicMatches().isEmpty()) {
                        // If it's size one, we can just compare the line
                        if (readAnnotation.getMergedTranscriptomicMatches().size() == 1) {
                            assertEquals(line, readAnnotation.output(), "Read annotation does not match the reference solution\n" +
                                    "Reference solution: " + line + "\n" +
                                    "Actual solution: " + readAnnotation.output() + LINE_SEPARATOR);
                        } else {

                            // ENSG00000215791,pseudogene:MERGED|ENSG00000160075,protein_coding:MERGED
                            checkMergedAndIntronicCorrectness(readAnnotation, parts, line, MERGED);
                        }
                    } else {
                        if (readAnnotation.getGenesThatInclude().size() == 1) {
                            assertEquals(line, readAnnotation.output(), "Read annotation does not match the reference solution\n" +
                                    "Reference solution: " + line + "\n" +
                                    "Actual solution: " + readAnnotation.output() + LINE_SEPARATOR);
                        } else {
                            // ENSG00000237094,lincRNA:INTRON|ENSG00000250575,pseudogene:INTRON
                            checkMergedAndIntronicCorrectness(readAnnotation, parts, line, INTRONIC);
                        }
                    }

                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        if (!readAnnotationMap.isEmpty()) {
            throw new IllegalStateException("There are read annotations that were not in the reference solution: " + readAnnotationMap.keySet());
        }
    }

    private static void checkMergedAndIntronicCorrectness(ReadAnnotation readAnnotation, String[] parts, String line, String TYPE) {
        String output = readAnnotation.output();
        String outputParts[] = output.split("\t");
        checkIfGeneralStatsMatch(outputParts, parts, line, output);
        String[] referenceGenes = parts[5].split("\\|");
        String[] actualGenes = outputParts[5].split("\\|");
        assertEquals(referenceGenes.length, actualGenes.length, "Number of genes does not match the reference solution\n" +
                "Reference solution: " + line + "\n" +
                "Actual solution: " + output + LINE_SEPARATOR);
        for (String refGen : referenceGenes) {
            String[] refGeneParts = refGen.split(":|,");
            String refGenNAME = refGeneParts[0];
            String refGenBioType = refGeneParts[1];
            boolean isGeneFound = false;
            for (String actualGene : actualGenes) {
                String[] actualGeneParts = actualGene.split(":|,");
                String actualGeneID = actualGeneParts[0];
                String actualGeneBioType = actualGeneParts[1];
                if (actualGeneID.equals(refGenNAME)) {
                    isGeneFound = true;
                    // BIotype
                    assertEquals(refGenBioType, actualGeneBioType, "Biotype does not match the reference solution\n" +
                            "Reference solution: " + line + "\n" +
                            "Actual solution: " + output + LINE_SEPARATOR);
                    // check if it says MERGED
                    assertEquals(TYPE, actualGeneParts[2], "merged does not match the reference solution\n" +
                            "Reference solution: " + line + "\n" +
                            "Actual solution: " + output + LINE_SEPARATOR);
                    break;
                }
            }
            if (!isGeneFound) {
                throw new IllegalStateException("Gene not found in the results:\n" +
                        "Reference solution: " + line + "\n" +
                        "Actual solution: " + output + LINE_SEPARATOR);
            }
        }
    }

    record ReferenceEntry(String bam, String gtf, Boolean strandness, String referenceSolution) {
    }

}
