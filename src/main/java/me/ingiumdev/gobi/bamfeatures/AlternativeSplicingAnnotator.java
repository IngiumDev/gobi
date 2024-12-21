package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import me.ingiumdev.gobi.gtf.ExonSkip;
import me.ingiumdev.gobi.gtf.GTFAnnotation;
import me.ingiumdev.gobi.gtf.structs.Gene;
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
        readsToAnnotate.stream().forEach(readPair -> {
            ReadAnnotation readAnnotation = processRead(readPair);
            if (readAnnotation != null) {
                // go through the gene matches etc
                for (Pair<Gene, List<Transcript>> match : readAnnotation.getTranscriptomicMatches()) {
                    String geneID = match.getFirst().getGeneID();
                    Set<String> matchTranscriptIds = match.getSecond().stream().map(Transcript::getTranscriptID).collect(Collectors.toSet());
                    List<ExonReadCount> exonReadCounts = readCounts.get(geneID);
                    if (exonReadCounts == null) {
                        continue;
                    }
                    for (ExonReadCount exonReadCount : exonReadCounts) {
                        boolean isWtMatched = exonReadCount.getWtTranscripts().equals(matchTranscriptIds);
                        boolean isSvMatched = exonReadCount.getSvTranscripts().equals(matchTranscriptIds);
                        if (isWtMatched) {
                            if (!isSvMatched) {
                                // Inclusion
                                if (readAnnotation.getAlignmentStart() <= exonReadCount.getSkippedExon().getEnd() && readAnnotation.getAlignmentEnd() >= exonReadCount.getSkippedExon().getStart()) {
                                    exonReadCount.incrementInclusionCount();
                                }
                            }
                        } else if (isSvMatched) {
                            // exclusion
                            if (readAnnotation.getAlignmentStart() < exonReadCount.getSkippedExon().getStart() && readAnnotation.getAlignmentEnd() > exonReadCount.getSkippedExon().getEnd()) {
                                exonReadCount.incrementExclusionCount();
                            }
                        }
                    }
                }
            }
        });
    }

    @Override
    protected ReadAnnotation processRead(SAMReadPair samReadPair) {
        SAMRecord first = samReadPair.getFirst();
        SAMRecord second = samReadPair.getSecond();
        ReadAnnotation readAnnotation = new ReadAnnotation(first.getReadName());

        readAnnotation.extractReadAlignmentStartEnd(first, second);

        List<Gene> resultGenes = forestManager.getGenesThatInclude(readAnnotation);
        if (!resultGenes.isEmpty()) {
            readAnnotation.extractReadIntervals(first, second);
            readAnnotation.calculateSplit();
            if (readAnnotation.areReadsConsistent()) {
                readAnnotation.setGenesThatInclude(resultGenes);
                if (readAnnotation.findTranscriptomicMatches()) {
                    return readAnnotation;
                }
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
