package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import me.ingiumdev.gobi.gtf.ExonSkip;
import me.ingiumdev.gobi.gtf.GTFAnnotation;
import me.ingiumdev.gobi.gtf.structs.Gene;
import me.ingiumdev.gobi.gtf.structs.Interval;
import me.ingiumdev.gobi.gtf.structs.Transcript;
import me.ingiumdev.gobi.gtf.treecollections.StrandUnspecificForest;
import me.ingiumdev.gobi.parsers.GTFParser;
import me.ingiumdev.gobi.readsimulator.Pair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static me.ingiumdev.gobi.bamfeatures.ExonReadCount.createReadCounts;

public class AlternativeSplicingAnnotator extends ReadAnnotator {
    List<ReadAnnotation> reads;
    Map<String, List<ExonReadCount>> readCounts;


    private AlternativeSplicingAnnotator(Builder builder) {
        samReader = builder.samReader;
        gtfFile = builder.gtfFile;
        outputFile = builder.outputFile;
        reads = new ArrayList<>();
    }


    @Override
    public void init() {
        forestManager = new StrandUnspecificForest();
        if (outputFile != null) {
            try {
                writer = new BufferedWriter(new FileWriter(outputFile));
                writer.write("gene" + "\t" + "exon" + "\t" + "num_incl_reads" + "\t" + "num_excl_reads" + "\t" + "num_total_reads" + "\t" + "psi");
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(String.valueOf(gtfFile));
        gtfAnnotation.getGenes().values().parallelStream().forEach(Gene::processIntrons);
        forestManager.init(gtfAnnotation);
        List<ExonSkip> exonSkips = ExonSkip.findExonSkippingEvents(gtfAnnotation);
        readCounts = createReadCounts(exonSkips);

    }

    @Override
    protected void processNewChromosome(String referenceName) {
        if (readsToAnnotate != null) {
            processReads();
        }
        readsToAnnotate = new ArrayList<>();
        lookup = new HashMap<>();
        currentChromosome = referenceName;
        forestManager.nextTree(referenceName);
    }

    private void processReads() {
        Set<String> seennames = new HashSet<>();

        Set<Integer> numbers = new HashSet<>(Set.of(66546, 66532, 66452, 66440, 66547, 66586, 66549, 66422, 66408, 66437, 66419, 66474, 66592, 66406, 66469, 66568, 66424, 66429, 66473, 66525, 66520));
        Set<Integer> found = new HashSet<>();
        readsToAnnotate.stream().forEach(readPair -> {
            ReadAnnotation readAnnotation = processRead(readPair);
            boolean seen = seennames.add(readPair.getFirst().getReadName());
            if (!seen) {
                System.out.println("Duplicate read name: " + readPair.getFirst().getReadName());
            }
            if (readAnnotation != null) {
                // go through the gene matches etc

                for (Pair<Gene, List<Transcript>> match : readAnnotation.getTranscriptomicMatches()) {
                    String geneID = match.getFirst().getGeneID();
                    Set<String> matchTranscriptIds = match.getSecond().stream()
                            .map(Transcript::getTranscriptID)
                            .collect(Collectors.toSet());
                    List<ExonReadCount> exonReadCounts = readCounts.get(geneID);
                    if (exonReadCounts == null) {
                        continue;
                    }
// TODO only print >= total above 0
                    for (ExonReadCount exonReadCount : exonReadCounts) {

                        boolean isWTmatched = exonReadCount.getWT_trans().equals(matchTranscriptIds);
                        boolean isSVmatched = exonReadCount.getSV_trans().equals(matchTranscriptIds);
                        if (exonReadCount.getSkippedExon().getStart() == 30702432 && exonReadCount.getSkippedExon().getEnd() == 30702470 && readAnnotation.getReadID().equals("77779")) {
//                            System.out.println("Inclusion: " + readAnnotation.getReadID() + " " + readBlock.getStart() + " " + readBlock.getEnd());
//                               found.add(Integer.valueOf(readAnnotation.getReadID()));
//                                                        System.out.println();
                        }
                        if (isWTmatched) {
                            if (!isSVmatched) {
                                // Inclusion
                                //check if there is an AB in the exon
                                int exonstart = exonReadCount.getSkippedExon().getStart();
                                int exonend = exonReadCount.getSkippedExon().getEnd();
                                if (false) {
                                    if (readAnnotation.getAlignmentStart() <= exonReadCount.getSkippedExon().getEnd() && readAnnotation.getAlignmentEnd() >= exonReadCount.getSkippedExon().getStart()) {
                                        exonReadCount.incrementExclusionCount();
                                    }
                                } else {
                                    for (Interval readBlock : readAnnotation.getCombinedRead()) {
                                        if (exonstart == 30702432 && exonend == 30702470 && readAnnotation.getReadID().equals("66592")) {
                                            System.out.println("Inclusion: " + readAnnotation.getReadID() + " " + readBlock.getStart() + " " + readBlock.getEnd());
                                            //   found.add(Integer.valueOf(readAnnotation.getReadID()));
                                        }
                                        if (readBlock.getStart() >= exonstart && readBlock.getEnd() <= exonend) {

                                            exonReadCount.incrementInclusionCount();
                                            if (exonstart == 30702432 && exonend == 30702470) {
//                                        System.out.println("Inclusion: " + readAnnotation.getReadID() + " " + readBlock.getStart() + " " + readBlock.getEnd());
                                                found.add(Integer.valueOf(readAnnotation.getReadID()));
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                        } else if (isSVmatched) {
                            if (readAnnotation.getAlignmentStart() < exonReadCount.getSkippedExon().getStart() && readAnnotation.getAlignmentEnd() > exonReadCount.getSkippedExon().getEnd()) {
                                exonReadCount.incrementExclusionCount();
                            }

                        }
                    }
                }
            }
        });

        if (numbers.removeAll(found)) {
            System.out.println("Found: " + found);
            System.out.println("Not found: " + numbers);
        }



    }

    @Override
    protected ReadAnnotation processRead(SAMReadPair samReadPair) {
        if (samReadPair.getFirst().getReadName().equals("77779")) {
            System.out.println();
        }
        var sr = samReadPair.getFirst();
        SAMRecord first = samReadPair.getFirst();
        SAMRecord second = samReadPair.getSecond();
        ReadAnnotation readAnnotation = new ReadAnnotation(first.getReadName());
        readAnnotation.extractReadIntervals(first, second);
        readAnnotation.extractReadAlignmentStartEnd(first, second);

        List<Gene> resultGenes = forestManager.getGenesThatInclude(readAnnotation);

        readAnnotation.calculateSplit();
        if (!resultGenes.isEmpty() && readAnnotation.areReadsConsistent()) {
            readAnnotation.setGenesThatInclude(resultGenes);
// remove split inconsistent check if not needed
            if (readAnnotation.findTranscriptomicMatches()) {
                return readAnnotation;
            }
        }
        return null;
    }

    @Override
    protected void processLastChromosome() {
        if (readsToAnnotate != null) {
            processReads();
        }
        // write
        readCounts.values().parallelStream().flatMap(Collection::stream).forEach(ExonReadCount::processPSI);
        readCounts.values().stream().flatMap(Collection::stream).forEach(exonReadCount -> {

        });
        if (outputFile != null) {


            for (var rc : readCounts.values()) {
                for (var erc : rc) {

                    try {

                        if ((erc.getExclusionCount() > 0 || erc.getInclusionCount() > 0)) {
                            writer.newLine();
                            writer.write(erc.output());
                        }

                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }


    }

    public static final class Builder {
        private SamReader samReader;
        private File gtfFile;
        private File outputFile;

        public Builder() {
        }

        public Builder setSamReader(SamReader samReader) {
            this.samReader = samReader;
            return this;
        }

        public Builder setGtfFile(File gtfFile) {
            this.gtfFile = gtfFile;
            return this;
        }

        public Builder setOutputFile(File outputFile) {
            this.outputFile = outputFile;
            return this;
        }

        public AlternativeSplicingAnnotator build() {
            return new AlternativeSplicingAnnotator(this);
        }
    }
}
