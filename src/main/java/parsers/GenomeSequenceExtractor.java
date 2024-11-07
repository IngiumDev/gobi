package parsers;

import gtf.structs.Interval;
import gtf.types.StrandDirection;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class GenomeSequenceExtractor {

    /* TESTCASES
    >ENST00000390522 cdna:pseudogene chromosome:GRCh37:14:22998580:22998641:1 gene:ENSG00000211874 gene_biotype:TR_J_pseudogene transcript_biotype:TR_J_pseudogene
CCAACCAGGCAGGGAACTGCTCTGATCTTTGGGAAGGGAACCCACCTTATCAGTGAGTTC
CA

>ENST00000504127 cdna:pseudogene chromosome:GRCh37:14:22956171:22956233:1 gene:ENSG00000248366 gene_biotype:TR_J_pseudogene transcript_biotype:TR_J_pseudogene
AGATGCGTGACAGCTATGAGAAGCTGATATTTGGAAAGGAGACATGACTAACTGTGAAGC
CAA


>ENST00000442641 cdna:pseudogene chromosome:GRCh37:14:22945296:22945352:1 gene:ENSG00000249446 gene_biotype:TR_J_pseudogene transcript_biotype:TR_J_pseudogene
TGAAGATCACCTAGATGCTCAACTTTGGGAAGGGGACTGAGTTAATTGTGAGCCTGG


>ENST00000509135 cdna:pseudogene chromosome:GRCh37:14:22950686:22950742:1 gene:ENSG00000250688 gene_biotype:TR_J_pseudogene transcript_biotype:TR_J_pseudogene
ACAAGTGCTGGTAATGCTCCTGTTGGGGAAAGGGGATGAGTACAAAAATAAATCCAA


>ENST00000415118 cdna:known chromosome:GRCh37:14:22907539:22907546:1 gene:ENSG00000223997 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
GAAATAGT

>ENST00000434970 cdna:known chromosome:GRCh37:14:22907999:22908007:1 gene:ENSG00000237235 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
CCTTCCTAC

>ENST00000448914 cdna:known chromosome:GRCh37:14:22918105:22918117:1 gene:ENSG00000228985 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
ACTGGGGGATACG
     */
    private final FileChannel fileChannel;
    private final Map<String, FastaIndexEntry> fastaIndex;

    public GenomeSequenceExtractor(String fastaFilePath, String faiFilePath) {
        try {
            this.fileChannel = FileChannel.open(Paths.get(fastaFilePath));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        this.fastaIndex = loadIndex(faiFilePath);
    }

    public static void main(String[] args) {
        GenomeSequenceExtractor genomeSequenceExtractor = new GenomeSequenceExtractor("C:\\Users\\Simon\\IdeaProjects\\gobi\\data\\readsimulator\\Homo_sapiens.GRCh37.75.dna.toplevel.fa", "C:\\Users\\Simon\\IdeaProjects\\gobi\\data\\readsimulator\\Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai");
        String sequence = genomeSequenceExtractor.getSequence("14", 106382685, 106382715, StrandDirection.REVERSE);
        System.out.println(sequence);
    }

    public static String reverseComplement(String input) {
        StringBuilder sb = new StringBuilder();
        for (int i = input.length() - 1; i >= 0; i--) {
            char c = input.charAt(i);
            switch (c) {
                case 'A' -> sb.append('T');
                case 'T' -> sb.append('A');
                case 'C' -> sb.append('G');
                case 'G' -> sb.append('C');
                default -> sb.append(c);
            }
        }
        return sb.toString();
    }

    private Map<String, FastaIndexEntry> loadIndex(String faiFilePath) {
        Map<String, FastaIndexEntry> index = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(faiFilePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                int fieldStart = 0;
                int tabCount = 0;

                String chrName = null;
                long byteOffset = 0;
                int lineBases = 0;
                int lineBytes = 0;
                int length = 0;

                for (int i = 0; i < line.length(); i++) {
                    if (line.charAt(i) == '\t') {
                        switch (tabCount) {
                            case 0 -> chrName = line.substring(fieldStart, i);
                            case 1 -> length = Integer.parseInt(line.substring(fieldStart, i));  // Sequence length
                            case 2 -> byteOffset = Long.parseLong(line.substring(fieldStart, i));
                            case 3 -> lineBases = Integer.parseInt(line.substring(fieldStart, i));
                        }
                        fieldStart = i + 1;
                        tabCount++;
                    }
                }
                if (fieldStart < line.length()) {
                    lineBytes = Integer.parseInt(line.substring(fieldStart));
                }
                index.put(chrName, new FastaIndexEntry(byteOffset, lineBases, lineBytes, length));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return index;
    }

    public String getSequence(String chr, int start, int end, StrandDirection strand) {
        String sequence = getSequenceBuilder(chr, start, end).toString();
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence);
        }
        return sequence;
    }

    public String getSequence(String chr, Interval interval, StrandDirection strand) {
        StringBuilder sequence = getSequenceBuilder(chr, interval);
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence.toString());
        }
        return sequence.toString();
    }

    public String getSplitSequence(String chr, TreeSet<Interval> intervals, StrandDirection strand) {
        StringBuilder sequence = getSplitSequenceBuilder(chr, intervals);
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence.toString());
        }
        return sequence.toString();
    }

    private String getSequence(String chr, int start, int end) {
        return getSequenceBuilder(chr, start, end).toString();
    }

    private StringBuilder getSequenceBuilder(String chr, Interval interval) {
        return getSequenceBuilder(chr, interval.getStart(), interval.getEnd());
    }


    private StringBuilder getSplitSequenceBuilder(String chr, TreeSet<Interval> intervals) {
        StringBuilder sequence = new StringBuilder();
        for (Interval interval : intervals) {
            sequence.append(getSequenceBuilder(chr, interval));
        }
        return sequence;
    }


    private StringBuilder getSequenceBuilder(String chr, int start, int end) {
        // Retrieve the index entry for the chromosome
        FastaIndexEntry entry = fastaIndex.get(chr);
        if (entry == null) {
            throw new IllegalArgumentException("Chromosome " + chr + " not found in FAI index.");
        }

        // Calculate the byte offset and length to map based on the start and end positions
        long byteOffset = entry.byteOffset() - 1 + (start / entry.lineBases()) * entry.lineBytes() + start % entry.lineBases();
        int sequenceLength = end - start + 1;
        try {
            MappedByteBuffer buffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, byteOffset, sequenceLength + (sequenceLength / entry.lineBases()) + entry.lineBytes());
            StringBuilder sequence = new StringBuilder();
            int lineBaseCount = 0;

            while (buffer.hasRemaining()) {
                byte b = buffer.get();
                if (b == '\n') {
                    // Skip the newline character
                    continue;
                }

                // Append the character to the sequence
                sequence.append((char) b);
                lineBaseCount++;

                // Stop if we've read enough sequence
                if (lineBaseCount >= sequenceLength) {
                    break;
                }
            }
            return sequence;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // https://manpages.ubuntu.com/manpages/xenial/man5/faidx.5.html#:~:text=An%20fai%20index%20file%20is,line%20LINEWIDTH%20The%20number%20of
    private record FastaIndexEntry(long byteOffset, int lineBases, int lineBytes, int length) {
    }

}
