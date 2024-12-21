package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import me.ingiumdev.gobi.gtf.ExonSkip;
import me.ingiumdev.gobi.gtf.GTFAnnotation;
import me.ingiumdev.gobi.gtf.structs.Gene;
import me.ingiumdev.gobi.gtf.treecollections.StrandUnspecificForest;
import me.ingiumdev.gobi.parsers.GTFParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class AlternativeSplicingAnnotator extends ReadAnnotator {
    List<ReadAnnotation> reads;

    private AlternativeSplicingAnnotator(Builder builder) {
        samReader = builder.samReader;
        gtfFile = builder.gtfFile;
        outputFile = builder.outputFile;
        reads = new ArrayList<>();
    }

//     TODO struct for the new exon skipping events


    @Override
    public void init() {
        forestManager = new StrandUnspecificForest();
        if (outputFile != null) {
            try {
                writer = new BufferedWriter(new FileWriter(outputFile));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(String.valueOf(gtfFile));
        gtfAnnotation.getGenes().values().parallelStream().forEach(Gene::processIntrons);
        List<ExonSkip> exonSkips = ExonSkip.findExonSkippingEvents(gtfAnnotation);
        forestManager.init(gtfAnnotation);
    }

    @Override
    protected void processNewChromosome(String referenceName) {
        if (readsToAnnotate != null) {
            // processing
        }
        readsToAnnotate = new ArrayList<>();
        lookup = new HashMap<>();
        currentChromosome = referenceName;
        forestManager.nextTree(referenceName);
    }

    @Override
    protected ReadAnnotation processRead(SAMReadPair samReadPair) {
        SAMRecord first = samReadPair.getFirst();
        SAMRecord second = samReadPair.getSecond();
        ReadAnnotation readAnnotation = new ReadAnnotation(first.getReadName());
        readAnnotation.extractReadAlignmentStartEnd(first, second);
        List<Gene> resultGenes = forestManager.getGenesThatInclude(readAnnotation);
        readAnnotation.extractReadIntervals(first, second);
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
