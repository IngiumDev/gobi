package readsimulator;

import gtf.GTFAnnotation;
import gtf.structs.Interval;
import gtf.structs.Transcript;
import parsers.GenomeSequenceExtractor;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static parsers.GTFParser.parseGTFForCounts;


public class ReadSimulator {
    private final int readLength;
    private final int meanFragmentLength;
    private final int fragmentLengthStandardDeviation;
    private final double mutationRate;
    private final GTFAnnotation gtfAnnotation;
    private final GenomeSequenceExtractor genomeSequenceExtractor;
    private final Map<String, Map<String, Integer>> readCounts;
    private final SplittableRandom random;
    private final String qualityScore;
    private final double L;


    public ReadSimulator(Builder builder) {
        this.readLength = builder.readLength;
        this.meanFragmentLength = builder.meanFragmentLength;
        this.fragmentLengthStandardDeviation = builder.fragmentLengthStandardDeviation;
        this.mutationRate = builder.mutationRate;
        this.gtfAnnotation = builder.gtfAnnotation;
        this.genomeSequenceExtractor = builder.genomeSequenceExtractor;
        this.readCounts = builder.readCounts;
        this.random = new SplittableRandom();
        this.qualityScore = "I".repeat(readLength);
        L = Math.exp(-(mutationRate * readLength));
    }

    private static void writeToReadMap(BufferedWriter mappingWriter, int readID, ReadPair rp) throws IOException {
        // Initialize StringBuilder for efficient string concatenation
        StringBuilder lineBuilder = new StringBuilder();

        // Append readID and tab
        lineBuilder.append(readID).append("\t");

        // Append sequence name, gene ID, and transcript ID
        lineBuilder.append(rp.getSeqName()).append("\t")
                .append(rp.getGeneID()).append("\t")
                .append(rp.getTranscriptID()).append("\t");

        // Append chromosomal coordinates for the first read
        TreeSet<Interval> firstCoordinates = rp.getFirst().getChromosomalCoordinates();
        appendCoordinates(lineBuilder, firstCoordinates);

        lineBuilder.append("\t");

        // Append chromosomal coordinates for the second read
        TreeSet<Interval> secondCoordinates = rp.getSecond().getChromosomalCoordinates();
        appendCoordinates(lineBuilder, secondCoordinates);

        lineBuilder.append("\t");

        // Append transcript coordinates for the first and second reads
        lineBuilder.append(rp.getFirst().getTranscriptCoordinates().getStart())
                .append("-")
                .append(rp.getFirst().getTranscriptCoordinates().getEnd() + 1)
                .append("\t");

        lineBuilder.append(rp.getSecond().getTranscriptCoordinates().getStart())
                .append("-")
                .append(rp.getSecond().getTranscriptCoordinates().getEnd() + 1)
                .append("\t");

        // Append mutated positions for the first read
        appendMutatedPositions(lineBuilder, rp.getFirst().getMutatedPositions());

        lineBuilder.append("\t");

        // Append mutated positions for the second read
        appendMutatedPositions(lineBuilder, rp.getSecond().getMutatedPositions());

        // Write the built line to the BufferedWriter
        mappingWriter.newLine();
        mappingWriter.write(lineBuilder.toString());
    }

    public int getReadLength() {
        return readLength;
    }

    public int getMeanFragmentLength() {
        return meanFragmentLength;
    }

    public int getFragmentLengthStandardDeviation() {
        return fragmentLengthStandardDeviation;
    }

    public double getMutationRate() {
        return mutationRate;
    }

    public GTFAnnotation getGtfAnnotation() {
        return gtfAnnotation;
    }

    public GenomeSequenceExtractor getGenomeSequenceExtractor() {
        return genomeSequenceExtractor;
    }

    public Map<String, Map<String, Integer>> getReadCounts() {
        return readCounts;
    }

    // Helper method to append coordinates
    private static void appendCoordinates(StringBuilder builder, TreeSet<Interval> coordinates) {
        boolean first = true;
        for (Interval interval : coordinates) {
            if (!first) {
                builder.append("|");
            }
            builder.append(interval.getStart()).append("-").append(interval.getEnd() + 1);
            first = false;
        }
    }



    public List<ReadPair> simulateReads() {
        return readCounts.entrySet().stream()
                .flatMap(geneEntry -> {
                    String geneId = geneEntry.getKey();
                    Map<String, Integer> transcriptMap = geneEntry.getValue();

//                    if (transcriptMap.size() > 1) {
//                        // Get the minimum start and maximum end positions for the gene
//                        int minStart = Integer.MAX_VALUE;
//                        int maxEnd = Integer.MIN_VALUE;
//                        for (Transcript transcript : gtfAnnotation.getGene(geneId).getTranscripts().values()) {
//                            minStart = Math.min(minStart, transcript.getExons().getFirst().getInterval().getStart());
//                            maxEnd = Math.max(maxEnd, transcript.getExons().getLast().getInterval().getEnd());
//                        }
//
//                        // readsimulator.Read the relevant section of the genome
//                        StringBuilder sequence = genomeSequenceExtractor.readSequence(gtfAnnotation.getGene(geneId).getSeqname(), minStart, maxEnd, gtfAnnotation.getGene(geneId).getStrand());
//
//                        // Simulate the reads
//                        int finalMinStart = minStart;
//                        return transcriptMap.entrySet().stream().flatMap(transcriptEntry -> {
//                            String transcriptId = transcriptEntry.getKey();
//                            int readCount = transcriptEntry.getValue();
//                            Transcript transcript = gtfAnnotation.getGene(geneId).getTranscript(transcriptId);
//                            // Basically, now we need to cut out the relevant section of the sequence for the transcript's exons
//                            int start = transcript.getExons().getFirst().getInterval().getStart();
//                            StringBuilder requiredSequence = new StringBuilder();
                    // TODO: fix this, we need to revcomp after cutting the sequence up, not before. Look into how to refactor this
//                            for (Exon exon : transcript.getExons()) {
//                                requiredSequence.append(sequence, exon.getInterval().getStart() - finalMinStart, exon.getInterval().getEnd() - finalMinStart + 1);
//                            }
//                            return simulateReadPairs(requiredSequence.toString(), transcript, readCount,transcript.getSeqname(),geneId , transcriptId).stream();
//
//                        });
//                    } else {
                    // Procedure for when there is only one transcript
                    return transcriptMap.entrySet().stream()
                            .flatMap(transcriptEntry -> {
                                String transcriptId = transcriptEntry.getKey();
                                int readCount = transcriptEntry.getValue();
                                Transcript transcript = gtfAnnotation.getGene(geneId).getTranscript(transcriptId);
                                String sequence = genomeSequenceExtractor.
                                        getSequenceForExonsInOneRead(gtfAnnotation.getGene(geneId).getSeqname(), transcript.getExons(), gtfAnnotation.getGene(geneId).getStrand());
                                return simulateReadPairs(sequence, transcript, readCount, transcript.getSeqname(), geneId, transcriptId).stream();
                            });
//                    }
                })
                .collect(Collectors.toList());
    }

    // Helper method to append mutated positions
    private static void appendMutatedPositions(StringBuilder builder, Set<Integer> positions) {
        boolean first = true;
        for (Integer position : positions) {
            if (!first) {
                builder.append(",");
            }
            builder.append(position);
            first = false;
        }
    }

    private void writeReadsToFASTQ(ReadPair rp, int readID, BufferedWriter fwWriter, BufferedWriter rwWriter) throws IOException {
        writeReads(readID, fwWriter, rp.getFirst());
        writeReads(readID, rwWriter, rp.getSecond());
    }

    /* readid	chr	gene	transcript	fw_regvec	rw_regvec	t_fw_regvec	t_rw_regvec	fw_mut	rw_mut
    *The simulator should output three files, two for the simulated paired-end sequences in
FASTQ format (fw.fastq, rw.fastq, one FASTQ file for the first read of a fragment, one
for the second; set the quality score to the maximum for all bases), and a tab separated file
(read.mappinginfo) with the following headers:
• readid: integer (starting from 0)
• chr id: chromosome
• gene id
• transcript id
• fw regvec: genomic region vector for the forward read (1-based end exclusive)
• rw regvec: genomic region vector for the reverse read (1-based end exclusive)
• t fw regvec: region vector for the forward read in transcript coordinates (0-based end
exclusive)
• t rw regvec: region vector for the reverse read in transcript coordinates (0-based end
exclusive)
• fw mut: mutated positions in the forward read (comma separated integer list)
• rw mut: mutated positions in the reverse read (comma separated integer list)
* The format for genomic region vectors is:
start1-end1(|startx-endx)+
Interval already has a toString method that outputs the interval but it's int the wrong format, it has a colon between the start and end but we need -. There can be multiple genome coordinates but only one transcript coordinate.
    *
    *
    * */
    public void writeReadCounts(String outputDir, List<ReadPair> readPairs) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir + "/read.mappinginfo"))) {
            writer.write("readid\tchr\tgene\ttranscript\tfw_regvec\trw_regvec\tt_fw_regvec\tt_rw_regvec\tfw_mut\trw_mut\n");
            int readId = 0;
            for (ReadPair readPair : readPairs) {
                Read first = readPair.getFirst();
                Read second = readPair.getSecond();
                writer.write(readId + "\t" + readPair.getSeqName() + "\t" + readPair.getGeneID() + "\t" + readPair.getTranscriptID() + "\t" +
                        first.getChromosomalCoordinates().stream().map(interval -> interval.getStart() + "-" + (interval.getEnd() + 1)).collect(Collectors.joining("|")) + "\t" +
                        second.getChromosomalCoordinates().stream().map(interval -> interval.getStart() + "-" + (interval.getEnd() + 1)).collect(Collectors.joining("|")) + "\t" +
                        first.getTranscriptCoordinates().getStart() + "-" + (first.getTranscriptCoordinates().getEnd() + 1) + "\t" +
                        second.getTranscriptCoordinates().getStart() + "-" + (second.getTranscriptCoordinates().getEnd() + 1) + "\t" +
                        first.getMutatedPositions().stream().map(String::valueOf).collect(Collectors.joining(",")) + "\t" +
                        second.getMutatedPositions().stream().map(String::valueOf).collect(Collectors.joining(","))
                );
                writer.newLine();
                readId++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeReads(String outputDir, List<ReadPair> readPairs) {
        // • fw.fastq
        //• rw.fastq
        //• one FASTQ file for the first read of a fragment, one for the second; set the quality score to the maximum for all bases
        try (BufferedWriter fwWriter = new BufferedWriter(new FileWriter(outputDir + "/fw.fastq"));
             BufferedWriter rwWriter = new BufferedWriter(new FileWriter(outputDir + "/rw.fastq"))) {
            int readId = 0;
            for (ReadPair readPair : readPairs) {
                Read first = readPair.getFirst();
                Read second = readPair.getSecond();
                writeReads(readId, fwWriter, first);

                writeReads(readId, rwWriter, second);
                readId++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void writeReads(int readID, BufferedWriter bw, Read read) throws IOException {
        //        bw.write("@" + readID);
//        bw.newLine();
//        bw.write(read.getSeq());
//        bw.newLine();
//        bw.write("+" + readID);
//        bw.newLine();
//        bw.write("I".repeat(read.getSeq().length()) + "\n");
        String stringBuilder = "@" + readID + "\n" +
                read.getSeq() + "\n" + "+" + readID + "\n" +
                qualityScore + "\n";
        //bw.newLine();
        bw.write(stringBuilder);
    }

    public void simulateAndWriteReads(String outputDir) {
//        long startTime;
//        long sequenceExtractionTime = 0;
//        long fragmentLengthSamplingTime = 0;
//        long readCreationTime = 0;
//        long mutationTime = 0;
//        long genomicPositionCalculationTime = 0;
//        long mappingInfoWritingTime = 0;
//        long fastqWritingTime = 0;
        try (BufferedWriter mappingWriter = new BufferedWriter(new FileWriter(outputDir + "/read.mappinginfo"));
             BufferedWriter fwWriter = new BufferedWriter(new FileWriter(outputDir + "/fw.fastq"));
             BufferedWriter rwWriter = new BufferedWriter(new FileWriter(outputDir + "/rw.fastq"))) {
            mappingWriter.write("readid\tchr\tgene\ttranscript\tfw_regvec\trw_regvec\tt_fw_regvec\tt_rw_regvec\tfw_mut\trw_mut");
            int readID = 0;
// TODO Move to byte[]
            for (String geneID : readCounts.keySet()) {
                for (String transcriptID : readCounts.get(geneID).keySet()) {
                    int readCount = readCounts.get(geneID).get(transcriptID);
                    Transcript transcript = gtfAnnotation.getGene(geneID).getTranscript(transcriptID);
//                    startTime = System.currentTimeMillis();
                    String sequence = genomeSequenceExtractor.getSequenceForExonsInOneRead(gtfAnnotation.getGene(geneID).getSeqname(), transcript.getExons(), gtfAnnotation.getGene(geneID).getStrand());
//                    sequenceExtractionTime += System.currentTimeMillis() - startTime;

                    for (int i = 0; i < readCount; i++) {
//                        startTime = System.currentTimeMillis();
                        int fragmentLength;
                        do {
                            fragmentLength = (int) Math.round(sampleFragmentLength());
                        } while (fragmentLength > sequence.length() || fragmentLength < readLength);
                        int diff = sequence.length() - fragmentLength;
                        int fragmentStart;
                        if (diff == 0) {
                            fragmentStart = 0;
                        } else {
                            fragmentStart = random.nextInt(diff);
                        }
//                        fragmentLengthSamplingTime += System.currentTimeMillis() - startTime;
//                        startTime = System.currentTimeMillis();
                        ReadPair rp = new ReadPair(sequence, fragmentStart, fragmentLength, readLength, transcript.getSeqname(), geneID, transcriptID, transcript.getStrand());
//                        readCreationTime += System.currentTimeMillis() - startTime;

//                        startTime = System.currentTimeMillis();

                        rp.mutateReadPairs(mutationRate, random, samplePoisson(), samplePoisson());
//                        mutationTime += System.currentTimeMillis() - startTime;

//                        startTime = System.currentTimeMillis();
                        rp.calculateGenomicPositions(transcript.getExons());
//                        genomicPositionCalculationTime += System.currentTimeMillis() - startTime;

//                        startTime = System.currentTimeMillis();
                        writeToReadMap(mappingWriter, readID, rp);
//                        mappingInfoWritingTime += System.currentTimeMillis() - startTime;
//                        startTime = System.currentTimeMillis();
                        writeReadsToFASTQ(rp, readID, fwWriter, rwWriter);
//                        fastqWritingTime += System.currentTimeMillis() - startTime;
                        readID++;
                    }
                }
            }
            // manually flush
//            startTime = System.currentTimeMillis();
            mappingWriter.flush();
//            mappingInfoWritingTime += System.currentTimeMillis() - startTime;
//            startTime = System.currentTimeMillis();
            fwWriter.flush();
            rwWriter.flush();
//            fastqWritingTime += System.currentTimeMillis() - startTime;
        } catch (IOException e) {
            e.printStackTrace();
        }
//        System.out.println("Sequence extraction time\t" + sequenceExtractionTime);
//        System.out.println("Fragment length sampling time\t" + fragmentLengthSamplingTime);
//        System.out.println("Read creation time\t" + readCreationTime);
//        System.out.println("Mutation time\t" + mutationTime);
//        System.out.println("Genomic position calculation time\t" + genomicPositionCalculationTime);
//        System.out.println("Mapping info writing time\t" + mappingInfoWritingTime);
//        System.out.println("Fastq writing time\t" + fastqWritingTime);
    }

    public double sampleFragmentLength() {
        return random.nextGaussian(meanFragmentLength, fragmentLengthStandardDeviation);
    }

    public int samplePoisson() {
        double p = 1.0;
        int k = 0;

        do {
            k++;
            p *= random.nextDouble();
        } while (p > L);

        return k - 1; // Return the number of events
    }

    public SplittableRandom getRandom() {
        return random;
    }

    private List<ReadPair> simulateReadPairs(String sequence, Transcript transcript, int readCount, String seqName, String geneID, String transcriptID) {

        List<ReadPair> readPairs = new ArrayList<>();
        for (int i = 0; i < readCount; i++) {
            int fragmentLength;
            do {
                fragmentLength = (int) Math.round(sampleFragmentLength());
                //  For simplicity, you may re-draw
                //the fragment length if the result is smaller than the read length or larger than the
                //transcript length.
            } while (fragmentLength > sequence.length() || fragmentLength < readLength);
            int fragmentStart = random.nextInt(sequence.length() - fragmentLength);
            ReadPair rp = new ReadPair(sequence, fragmentStart, fragmentLength, readLength, seqName, geneID, transcriptID, transcript.getStrand());
            // rp.mutateReadPairs(mutationRate, random, readLength);
            rp.calculateGenomicPositions(transcript.getExons());
            readPairs.add(rp);
        }
        return readPairs;
    }

    public static class Builder {
        private int readLength;
        private int meanFragmentLength;
        private int fragmentLengthStandardDeviation;
        private double mutationRate;
        private GTFAnnotation gtfAnnotation;
        private GenomeSequenceExtractor genomeSequenceExtractor;
        private Map<String, Map<String, Integer>> readCounts;

        public Builder setReadLength(int readLength) {
            this.readLength = readLength;
            return this;
        }

        public Builder setMeanFragmentLength(int meanFragmentLength) {
            this.meanFragmentLength = meanFragmentLength;
            return this;
        }

        public Builder setFragmentLengthStandardDeviation(int fragmentLengthStandardDeviation) {
            this.fragmentLengthStandardDeviation = fragmentLengthStandardDeviation;
            return this;
        }

        public Builder setMutationRate(double mutationRate) {
            this.mutationRate = mutationRate;
            return this;
        }

        public Builder setGtfAnnotation(String gtfPath) {
            this.gtfAnnotation = parseGTFForCounts(gtfPath, readCounts);
//            this.gtfAnnotation = parseGTF(gtfPath);
            return this;
        }

        public Builder setGenomeSequenceExtractor(GenomeSequenceExtractor genomeSequenceExtractor) {
            this.genomeSequenceExtractor = genomeSequenceExtractor;
            return this;
        }

        public Builder setReadCounts(Map<String, Map<String, Integer>> readCounts) {
            this.readCounts = readCounts;
            return this;
        }

        public ReadSimulator build() {
            return new ReadSimulator(this);
        }

    }
}
