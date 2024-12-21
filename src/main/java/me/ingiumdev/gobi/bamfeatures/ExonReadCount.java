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
    private final Set<String> svTranscripts;
    private final Set<String> wtTranscripts;
    private int totalCount;
    private double psi;

    private ExonReadCount(Builder builder) {
        geneID = builder.geneID;
        skippedExon = builder.skippedExon;
        intron = builder.intron;
        svTranscripts = builder.svTranscripts;
        wtTranscripts = builder.wtTranscripts;
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
                        .setSvTranscripts(SV_trans)
                        .setWtTranscripts(WT_trans)
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
        return geneID + "\t" + skippedExon.toStringDash() + "\t" + inclusionCount + "\t" + exclusionCount + "\t" + totalCount + "\t" + psi;
    }

    public String getGeneID() {
        return geneID;
    }

    public Interval getSkippedExon() {
        return skippedExon;
    }

    public Interval getIntron() {
        return intron;
    }

    public int getInclusionCount() {
        return inclusionCount;
    }

    public int getExclusionCount() {
        return exclusionCount;
    }

    public Set<String> getSvTranscripts() {
        return svTranscripts;
    }

    public Set<String> getWtTranscripts() {
        return wtTranscripts;
    }

    public int getTotalCount() {
        return totalCount;
    }

    public double getPsi() {
        return psi;
    }

    public static final class Builder {
        private String geneID;
        private Interval skippedExon;
        private Interval intron;
        private Set<String> svTranscripts;
        private Set<String> wtTranscripts;

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

        public Builder setSvTranscripts(Set<String> svTranscripts) {
            this.svTranscripts = svTranscripts;
            return this;
        }

        public Builder setWtTranscripts(Set<String> wtTranscripts) {
            this.wtTranscripts = wtTranscripts;
            return this;
        }

        public ExonReadCount build() {
            return new ExonReadCount(this);
        }
    }
}
