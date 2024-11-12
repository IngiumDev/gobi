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
    // TODO: Check if genomic coordinates match by when we generate a read and get the genomic coordinates, we can just get the read from these coordinates and check if it's the same, and we know our getsequence is correct
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
        String sequence = readSequence(chr, start, end).toString();
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence);
        }
        return sequence;
    }

    public String getSequence(String chr, Interval interval, StrandDirection strand) {
        StringBuilder sequence = readSequence(chr, interval);
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence.toString());
        }
        return sequence.toString();
    }

    public String getSequenceForIntervals(String chr, TreeSet<Interval> intervals, StrandDirection strand) {
        StringBuilder sequence = readSequenceForIntervals(chr, intervals);
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence.toString());
        }
        return sequence.toString();
    }

    // TODO: Get sequence for exons
    public String getSequenceForIntervalsInOneRead(String chr, TreeSet<Interval> intervals, StrandDirection strand) {
        StringBuilder sequence = readSequenceForIntervalsInOneRead(chr, intervals);
        if (strand == StrandDirection.REVERSE) {
            return reverseComplement(sequence.toString());
        }
        return sequence.toString();
    }

    private StringBuilder readSequenceForIntervalsInOneRead(String chr, TreeSet<Interval> intervals) {
        // Get a read that spans all intervals, then cut out the non required parts
        int start = intervals.first().getStart();
        int end = intervals.last().getEnd();
        StringBuilder sequence = readSequence(chr, start, end);
        StringBuilder result = new StringBuilder();
        for (Interval interval : intervals) {
            result.append(sequence, interval.getStart() - start, interval.getEnd() - start + 1);
        }
        return result;
    }

    private String getSequence(String chr, int start, int end) {
        return readSequence(chr, start, end).toString();
    }

    private StringBuilder readSequence(String chr, Interval interval) {
        return readSequence(chr, interval.getStart(), interval.getEnd());
    }


    private StringBuilder readSequenceForIntervals(String chr, TreeSet<Interval> intervals) {
        StringBuilder sequence = new StringBuilder();
        for (Interval interval : intervals) {
            sequence.append(readSequence(chr, interval));
        }
        return sequence;
    }

    public StringBuilder readSequence(String chr, int start, int end, StrandDirection direction) {
        StringBuilder sequence = readSequence(chr, start, end);
        if (direction == StrandDirection.REVERSE) {
            // TODO Revcomplement with stringbuilder
            return new StringBuilder(reverseComplement(sequence.toString()));
        }
        return sequence;
    }

    // ENSG00000164465	ENST00000368503	112500 starts at the last character in the line
    public StringBuilder readSequence(String chr, int start, int end) {
        // Retrieve the index entry for the chromosome
        FastaIndexEntry entry = fastaIndex.get(chr);
        if (entry == null) {
            throw new IllegalArgumentException("Chromosome " + chr + " not found in FAI index.");
        }

        // Calculate the byte offset and length to map based on the start and end positions
        long byteOffset = entry.byteOffset() + (long) ((start - 1) / entry.lineBases()) * entry.lineBytes() + ((start - 1) % entry.lineBases());
        int sequenceLength = end - start + 1;
        try {
            MappedByteBuffer buffer = fileChannel.map(FileChannel.MapMode.READ_ONLY,
                    byteOffset,
                    sequenceLength + (sequenceLength / entry.lineBases()) + entry.lineBytes());
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

    public boolean doesChromosomeExist(String chrToCheck) {
        return fastaIndex.containsKey(chrToCheck);
    }

    // https://manpages.ubuntu.com/manpages/xenial/man5/faidx.5.html#:~:text=An%20fai%20index%20file%20is,line%20LINEWIDTH%20The%20number%20of
    private record FastaIndexEntry(long byteOffset, int lineBases, int lineBytes, int length) {
    }
}
