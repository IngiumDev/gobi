package me.ingiumdev.gobi.bamfeatures;

import me.ingiumdev.gobi.gtf.structs.Interval;

import java.util.Set;

public class ExonReadCount {
    private final String geneID;
    private final Interval skippedExon;
    private final Interval intron;
    private final int inclusionCount;
    private final int exclusionCount;
    private int totalCount;
    private double psi;
    private final Set<String> SV_trans;
    private final Set<String> WT_trans;

    private ExonReadCount(Builder builder) {
        geneID = builder.geneID;
        skippedExon = builder.skippedExon;
        intron = builder.intron;
        SV_trans = builder.SV_trans;
        WT_trans = builder.WT_trans;
        inclusionCount = 0;
        exclusionCount = 0;
    }

    public void processPSI() {
        totalCount = inclusionCount + exclusionCount;
        psi = (double) inclusionCount / totalCount;
    }


    public static final class Builder {
        private String geneID;
        private Interval skippedExon;
        private Interval intron;
        private Set<String> SV_trans;
        private Set<String> WT_trans;

        public Builder() {
        }


        public Builder setGeneID(String geneID) {
            this.geneID = geneID;
            return this;
        }

        public Builder setSkippedExon(Interval skippedExon) {
            this.skippedExon = skippedExon;
            return this;
        }

        public Builder setIntron(Interval intron) {
            this.intron = intron;
            return this;
        }

        public Builder setSV_trans(Set<String> SV_trans) {
            this.SV_trans = SV_trans;
            return this;
        }

        public Builder setWT_trans(Set<String> WT_trans) {
            this.WT_trans = WT_trans;
            return this;
        }

        public ExonReadCount build() {
            return new ExonReadCount(this);
        }
    }
}
