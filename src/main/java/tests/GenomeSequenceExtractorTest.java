package tests;

import gtf.GTFAnnotation;
import org.junit.jupiter.api.Test;
import parsers.GTFParser;
import parsers.GenomeSequenceExtractor;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

public class GenomeSequenceExtractorTest {


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