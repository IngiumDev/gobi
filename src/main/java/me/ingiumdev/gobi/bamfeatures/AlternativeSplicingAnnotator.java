package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SamReader;

import java.io.File;

public class AlternativeSplicingAnnotator extends ReadAnnotator {
    private AlternativeSplicingAnnotator(Builder builder) {
        samReader = builder.samReader;
        gtfFile = builder.gtfFile;
        outputFile = builder.outputFile;
    }

//     TODO struct for the new exon skipping events

    @Override
    public void annotateReads() {

    }

    @Override
    public void init() {

    }

    @Override
    protected void processNewChromosome(String referenceName) {

    }

    @Override
    protected ReadAnnotation processRead(SAMReadPair samReadPair) {
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
