package me.ingiumdev.gobi.bamfeatures;

import me.ingiumdev.gobi.gtf.ExonSkip;
import me.ingiumdev.gobi.gtf.structs.Interval;

import java.util.*;

public class ExonReadCount {
    private final String geneID;
    private final Interval skippedExon;
    private final Interval intron;
    private int inclusionCount;
    private int exclusionCount;
    private final Set<String> SV_trans;
    private final Set<String> WT_trans;
    private int totalCount;
    private double psi;

    private ExonReadCount(Builder builder) {
        geneID = builder.geneID;
        skippedExon = builder.skippedExon;
        intron = builder.intron;
        SV_trans = builder.SV_trans;
        WT_trans = builder.WT_trans;
        inclusionCount = 0;
        exclusionCount = 0;
    }

    public static Map<String, List<ExonReadCount>> createReadCounts(List<ExonSkip> exonSkips) {
        Map<String, List<ExonReadCount>> readCounts = new HashMap<>();

        for (ExonSkip exonSkip : exonSkips) {
            String geneID = exonSkip.getId();
            if (!readCounts.containsKey(geneID)) {
                readCounts.put(geneID, new ArrayList<>());
            }
            List<ExonReadCount> exonReadCounts = readCounts.get(geneID);
            Interval intron = exonSkip.getSV();
            Set<String> SV_trans = exonSkip.getSV_trans();
            Set<String> WT_trans = exonSkip.getWT_trans();
            for (Interval exon : exonSkip.getExonSkipped()) {
                ExonReadCount exonReadCount = new ExonReadCount.Builder()
                        .setGeneID(geneID)
                        .setSkippedExon(exon)
                        .setIntron(intron)
                        .setSV_trans(SV_trans)
                        .setWT_trans(WT_trans)
                        .build();
                exonReadCounts.add(exonReadCount);
            }
        }
        return readCounts;
    }

    public void processPSI() {
        totalCount = inclusionCount + exclusionCount;
        psi = (double) inclusionCount / totalCount;
    }

    public synchronized void incrementInclusionCount() {
        inclusionCount++;
    }

    public synchronized void incrementExclusionCount() {
        exclusionCount++;
    }

    public String output() {
        return geneID + "\t" + skippedExon + "\t" + inclusionCount + "\t" + exclusionCount + "\t" + totalCount + "\t" + psi;
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
