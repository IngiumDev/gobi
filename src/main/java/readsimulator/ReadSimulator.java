package readsimulator;

import gtf.GTFAnnotation;
import gtf.structs.Exon;
import gtf.structs.Transcript;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.NormalDistribution;
import parsers.GenomeSequenceExtractor;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;


public class ReadSimulator {
    NormalDistribution normalDistribution;
    private final int readLength;
    private final int meanFragmentLength;
    private final int fragmentLengthStandardDeviation;
    private final double mutationRate;
    private final GTFAnnotation gtfAnnotation;
    private final GenomeSequenceExtractor genomeSequenceExtractor;
    private final Map<String, Map<String, Integer>> readCounts;
    private final ContinuousDistribution.Sampler sampler;
    private final Random random;
    private final UniformRandomProvider rng;

    public ReadSimulator(Builder builder) {
        this.readLength = builder.readLength;
        this.meanFragmentLength = builder.meanFragmentLength;
        this.fragmentLengthStandardDeviation = builder.fragmentLengthStandardDeviation;
        this.mutationRate = builder.mutationRate;
        this.gtfAnnotation = builder.gtfAnnotation;
        this.genomeSequenceExtractor = builder.genomeSequenceExtractor;
        this.readCounts = builder.readCounts;

        this.random = new Random();
        this.normalDistribution = NormalDistribution.of(meanFragmentLength, fragmentLengthStandardDeviation);
        this.rng = RandomSource.create(RandomSource.MT);
        this.sampler = normalDistribution.createSampler(rng);
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
                                String sequence = genomeSequenceExtractor.getSequenceForExonsInOneRead(gtfAnnotation.getGene(geneId).getSeqname(), transcript.getExons(), gtfAnnotation.getGene(geneId).getStrand());
                                return simulateReadPairs(sequence, transcript, readCount, transcript.getSeqname(), geneId, transcriptId).stream();
                            });
//                    }
                })
                .collect(Collectors.toList());
    }

    private List<ReadPair> simulateReadPairs(String sequence, Transcript transcript, int readCount, String seqName, String geneID, String transcriptID) {

        List<ReadPair> readPairs = new ArrayList<>();
        for (int i = 0; i < readCount; i++) {
            int fragmentLength;
            do {
                fragmentLength = (int) Math.round(sampler.sample());
                //  For simplicity, you may re-draw
                //the fragment length if the result is smaller than the read length or larger than the
                //transcript length.
            } while (fragmentLength > sequence.length());
            fragmentLength=Math.max(fragmentLength,readLength);
            int fragmentStart = random.nextInt(sequence.length() - fragmentLength);
            ReadPair rp = new ReadPair(sequence, fragmentStart, fragmentLength, readLength, seqName, geneID, transcriptID, transcript.getStrand());
            rp.mutateReadPairs(mutationRate, random);
            rp.calculateGenomicPositions(transcript.getExons().stream()
                    .map(Exon::getInterval)
                    .collect(Collectors.toCollection(TreeSet::new)));
            readPairs.add(rp);
        }
        return readPairs;
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
                        first.getChromosomalCoordinates().stream().map(interval -> interval.getStart() + "-" + (interval.getEnd()+1)).collect(Collectors.joining("|")) + "\t" +
                        second.getChromosomalCoordinates().stream().map(interval -> interval.getStart() + "-" + (interval.getEnd()+1)).collect(Collectors.joining("|")) + "\t" +
                        first.getTranscriptCoordinates().getStart() + "-" + (first.getTranscriptCoordinates().getEnd()+1) + "\t" +
                        second.getTranscriptCoordinates().getStart() + "-" + (second.getTranscriptCoordinates().getEnd()+1) + "\t" +
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
                fwWriter.write("@" + readId);
                fwWriter.newLine();
                fwWriter.write(first.getSeq());
                fwWriter.newLine();
                fwWriter.write("+" + readId);
                fwWriter.newLine();
                fwWriter.write("I".repeat(first.getSeq().length()));
                fwWriter.newLine();

                rwWriter.write("@" + readId);
                rwWriter.newLine();
                rwWriter.write(second.getSeq());
                rwWriter.newLine();
                rwWriter.write("+" + readId);
                rwWriter.newLine();
                rwWriter.write("I".repeat(second.getSeq().length()));
                rwWriter.newLine();
                readId++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
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

        public Builder setGtfAnnotation(GTFAnnotation gtfAnnotation) {
            this.gtfAnnotation = gtfAnnotation;
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
