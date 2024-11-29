package tests;

import gtf.GTFAnnotation;
import gtf.structs.Transcript;
import org.junit.jupiter.api.Test;
import parsers.GTFParser;
import parsers.GenomeSequenceExtractor;
import readsimulator.ReadSimulator;
import readsimulator.SimulatedReadPair;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;
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


        ReadSimulator readSimulator = new ReadSimulator.Builder()
                .setReadLength(75)
                .setMeanFragmentLength(200)
                .setFragmentLengthStandardDeviation(80)
                .setMutationRate(1.0 / 100)
                .setReadCounts(readCountsFile(basePath.resolve("readcounts.simulation").toString()))
                .setGtfAnnotation(basePath.resolve("Homo_sapiens.GRCh37.75.gtf").toString())
                .setGenomeSequenceExtractor(genomeSequenceExtractor)
                .build();
        for (String geneID : readSimulator.getReadCounts().keySet()) {
            for (String transcriptID : readSimulator.getReadCounts().get(geneID).keySet()) {
                int readCount = readSimulator.getReadCounts().get(geneID).get(transcriptID);
                Transcript transcript = readSimulator.getGtfAnnotation().getGene(geneID).getTranscript(transcriptID);
                String sequence = genomeSequenceExtractor.getSequenceForExonsInOneRead(
                        readSimulator.getGtfAnnotation().getGene(geneID).getSeqname(),
                        transcript.getExons(),
                        readSimulator.getGtfAnnotation().getGene(geneID).getStrand()
                );
                String reverseComplement = reverseComplement(sequence);
                for (int i = 0; i < readCount; i++) {
                    int fragmentLength;
                    do {
                        fragmentLength = (int) Math.round(readSimulator.sampleFragmentLength());
                    } while (fragmentLength > sequence.length() || fragmentLength < readSimulator.getReadLength());

                    int diff = sequence.length() - fragmentLength;
                    int fragmentStart = (diff == 0) ? 0 : readSimulator.getRandom().nextInt(diff);

                    SimulatedReadPair rp = new SimulatedReadPair(sequence, fragmentStart, fragmentLength,
                            readSimulator.getReadLength(), transcript.getSeqname(),
                            geneID, transcriptID, transcript.getStrand(), reverseComplement);

                    rp.mutateReadPairs(readSimulator.getMutationRate(), readSimulator.getRandom(), readSimulator.samplePoisson(), readSimulator.samplePoisson());
                    rp.calculateGenomicPositions(transcript.getExons());

                    // Get the original and mutated sequences
                    String originalSequenceFW = genomeSequenceExtractor.getSequenceForIntervalsInOneRead(
                            transcript.getSeqname(), rp.getFirst().getChromosomalCoordinates(), transcript.getStrand()
                    );
                    String originalSequenceRV = reverseComplement(genomeSequenceExtractor.getSequenceForIntervalsInOneRead(
                            transcript.getSeqname(), rp.getSecond().getChromosomalCoordinates(), transcript.getStrand()
                    ));

                    String mutatedSequenceFW = rp.getFirst().getSeq();
                    String mutatedSequenceRV = rp.getSecond().getSeq();

                    // Convert mutated positions to sets for faster lookup
                    Set<Integer> mutatedPositionsFW = new HashSet<>(rp.getFirst().getMutatedPositions());
                    Set<Integer> mutatedPositionsRV = new HashSet<>(rp.getSecond().getMutatedPositions());

                    // Check mutations and unchanged sequences in a single loop
                    int lengthFW = mutatedSequenceFW.length();
                    int lengthRV = mutatedSequenceRV.length();
                    int maxLength = Math.max(lengthFW, lengthRV);

                    for (int j = 0; j < maxLength; j++) {
                        // Check forward read mutations and unchanged positions
                        if (j < lengthFW) {
                            if (mutatedPositionsFW.contains(j)) {
                                assertNotEquals(mutatedSequenceFW.charAt(j), originalSequenceFW.charAt(j),
                                        "Mutation not found at position " + j + " in forward read");
                            } else {
                                assertEquals(mutatedSequenceFW.charAt(j), originalSequenceFW.charAt(j),
                                        "Mismatch at position " + j + " in forward read");
                            }
                        }

                        // Check reverse read mutations and unchanged positions
                        if (j < lengthRV) {
                            if (mutatedPositionsRV.contains(j)) {
                                assertNotEquals(mutatedSequenceRV.charAt(j), originalSequenceRV.charAt(j),
                                        "Mutation not found at position " + j + " in reverse read");
                            } else {
                                assertEquals(mutatedSequenceRV.charAt(j), originalSequenceRV.charAt(j),
                                        "Mismatch at position " + j + " in reverse read");
                            }
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