package tests;

import gtf.GTFAnnotation;
import gtf.structs.Transcript;
import org.junit.jupiter.api.Test;
import parsers.GTFParser;
import parsers.GenomeSequenceExtractor;
import readsimulator.ReadPair;
import readsimulator.ReadSimulator;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.*;
import static parsers.GenomeSequenceExtractor.reverseComplement;
import static runners.ReadSimulatorRunner.readCountsFile;

public class ReadSimulatorTest {


    @Test
    public void testGenomicCoordinates() {
        // Test if the genomic coordinates are calculated correctly for a read
        Path basePath = Paths.get("data", "readsimulator");

        GenomeSequenceExtractor genomeSequenceExtractor = new GenomeSequenceExtractor(
                basePath.resolve("Homo_sapiens.GRCh37.75.dna.toplevel.fa").toString(),
                basePath.resolve("Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai").toString()
        );
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(
                basePath.resolve("Homo_sapiens.GRCh37.75.gtf").toString()
        );

        ReadSimulator readSimulator = new ReadSimulator.Builder()
                .setReadLength(75)
                .setMeanFragmentLength(200)
                .setFragmentLengthStandardDeviation(80)
                .setMutationRate(1.0 / 100)
                .setGtfAnnotation(gtfAnnotation)
                .setGenomeSequenceExtractor(genomeSequenceExtractor)
                .setReadCounts(readCountsFile(basePath.resolve("readcounts.simulation").toString()))
                .build();
        for (String geneID : readSimulator.getReadCounts().keySet()) {
            for (String transcriptID : readSimulator.getReadCounts().get(geneID).keySet()) {
                int readCount = readSimulator.getReadCounts().get(geneID).get(transcriptID);
                Transcript transcript = gtfAnnotation.getGene(geneID).getTranscript(transcriptID);
                String sequence = genomeSequenceExtractor.getSequenceForExonsInOneRead(gtfAnnotation.getGene(geneID).getSeqname(), transcript.getExons(), gtfAnnotation.getGene(geneID).getStrand());
                for (int i = 0; i < readCount; i++) {
                    int fragmentLength;
                    do {
                        fragmentLength = (int) Math.round(readSimulator.sampleFragmentLength());
                    } while (fragmentLength > sequence.length());
                    fragmentLength = Math.max(fragmentLength, readSimulator.getReadLength());
                    int diff = sequence.length() - fragmentLength;
                    int fragmentStart;
                    if (diff == 0) {
                        fragmentStart = 0;
                    } else {
                        fragmentStart = readSimulator.getRandom().nextInt(diff);
                    }
                    ReadPair rp = new ReadPair(sequence, fragmentStart, fragmentLength, readSimulator.getReadLength(), transcript.getSeqname(), geneID, transcriptID, transcript.getStrand());
                    rp.mutateReadPairs(readSimulator.getMutationRate(), readSimulator.getRandom(), readSimulator.getRng());
                    rp.calculateGenomicPositions(transcript.getExons());

                    // Check if the genomic coordinates are calculated correctly, get the original sequence and compare
                    String originalSequenceFW = genomeSequenceExtractor.getSequenceForIntervalsInOneRead(transcript.getSeqname(), rp.getFirst().getChromosomalCoordinates(), transcript.getStrand());
                    String originalSequenceRV = reverseComplement(genomeSequenceExtractor.getSequenceForIntervalsInOneRead(transcript.getSeqname(), rp.getSecond().getChromosomalCoordinates(), transcript.getStrand()));

                    // We check if the mutations are actually present in the sequence, and the rest of the sequence is the same
                    String mutatedSequenceFW = rp.getFirst().getSeq();
                    String mutatedSequenceRV = rp.getSecond().getSeq();

                    for (int pos : rp.getFirst().getMutatedPositions()) {
                        assertNotEquals(mutatedSequenceFW.charAt(pos), originalSequenceFW.charAt(pos), "Mutation not found at position " + pos + " in forward read");
                    }

                    for (int pos : rp.getSecond().getMutatedPositions()) {
                        assertNotEquals(mutatedSequenceRV.charAt(pos), originalSequenceRV.charAt(pos), "Mutation not found at position " + pos + " in reverse read");
                    }

                    // Check if the rest of the sequence is the same
                    for (int j = 0; j < mutatedSequenceFW.length(); j++) {
                        if (!rp.getFirst().getMutatedPositions().contains(j)) {
                            assertEquals(mutatedSequenceFW.charAt(j), originalSequenceFW.charAt(j), "Mismatch at position " + j + " in forward read");
                        }
                    }

                    for (int j = 0; j < mutatedSequenceRV.length(); j++) {
                        if (!rp.getSecond().getMutatedPositions().contains(j)) {
                            assertEquals(mutatedSequenceRV.charAt(j), originalSequenceRV.charAt(j), "Mismatch at position " + j + " in reverse read");
                        }
                    }

                }
            }
        }


    }

    @Test
    public void testGenomeSequenceExtractor() {
        // Use relative paths from the project root
        Path basePath = Paths.get("data", "readsimulator");

        GenomeSequenceExtractor genomeSequenceExtractor = new GenomeSequenceExtractor(
                basePath.resolve("Homo_sapiens.GRCh37.75.dna.toplevel.fa").toString(),
                basePath.resolve("Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai").toString()
        );
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(
                basePath.resolve("Homo_sapiens.GRCh37.75.gtf").toString()
        );
        String transcriptomePath = basePath.resolve("Homo_sapiens.GRCh37.75.cdna.all.fa.gz").toString();

        try (FileInputStream fileInputStream = new FileInputStream(transcriptomePath);
             GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
             BufferedReader br = new BufferedReader(new InputStreamReader(gzipInputStream))) {

            String line;
            String id = null;
            StringBuilder seq = new StringBuilder();
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (!seq.isEmpty()) {
                        // Check the sequence
                        checkIfReferenceSequenceMatches(id, seq.toString(), genomeSequenceExtractor, gtfAnnotation);
                        seq.setLength(0);
                    }
                    id = line.substring(1);
                } else {
                    seq.append(line);
                }
            }
            checkIfReferenceSequenceMatches(id,
                    seq.toString(),
                    genomeSequenceExtractor,
                    gtfAnnotation);
        } catch (IOException e) {
            fail("Test failed due to IOException: " + e.getMessage());
        }
    }

    private void checkIfReferenceSequenceMatches(String id, String referenceSequence, GenomeSequenceExtractor genomeSequenceExtractor, GTFAnnotation gtfAnnotation) {
        String[] parts = id.split(" ");
        String transcriptID = parts[0];
        String geneID = parts[3].split(":")[1];

        if (!genomeSequenceExtractor.doesChromosomeExist(gtfAnnotation.getGene(geneID).getSeqname())) {
            return;
        }

        String sequence = genomeSequenceExtractor.getSequenceForExonsInOneRead(
                gtfAnnotation.getGene(geneID).getSeqname(),
                gtfAnnotation.getGene(geneID).getTranscript(transcriptID).getExons(),
                gtfAnnotation.getGene(geneID).getStrand()
        );

        // Use assertions to validate the sequences
        assertEquals(referenceSequence, sequence, "Sequences do not match for transcript: " + transcriptID);
    }
}