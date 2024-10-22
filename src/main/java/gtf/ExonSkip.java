package gtf;

import gtf.structs.Exon;
import gtf.structs.Gene;
import gtf.structs.Interval;
import gtf.structs.Transcript;
import gtf.types.StrandDirection;
import parsers.GTFParser;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

public class ExonSkip {

    //gene_id
    private String id;
    // gene symbol
    private String symbol;
    // chromosome
    private String chr;
    // strand
    private StrandDirection strand;
    // number of annotated CDS in the gene
    private int nprots;
    // number of annotated transcripts in the gene
    private int ntrans;
    // the SV intron as start:end
    private Interval SV;
    // all the WT introns within the SV intron separated by | as start:end
    private TreeSet<Interval> WT;
    // ids of the SV CDS, separated by | (protein_id)
    private Set<String> SV_prots;
    // ids of the WT CDS, separated by | (protein_id)
    private Set<String> WT_prots;
    // the minimal number of skipped exons in any WT/SV pair
    private int min_skipped_exon;
    // the maximum number of skipped exons in any WT/SV pair
    private int max_skipped_exon;
    // the minimal number of skipped bases (joint length of skipped exons) in any WT/SV pair
    private int min_skipped_bases;
    // the maximum number of skipped bases (joint length of skipped exons) in any WT/SV pair
    private int max_skipped_bases;

    private ExonSkip(Builder builder) {
        this.id = builder.id;
        this.symbol = builder.symbol;
        this.chr = builder.chr;
        this.strand = builder.strand;
        this.nprots = builder.nprots;
        this.ntrans = builder.ntrans;
        this.SV = builder.SV;
        this.WT = builder.WT;
        this.SV_prots = builder.SV_prots;
        this.WT_prots = builder.WT_prots;
        this.min_skipped_exon = builder.min_skipped_exon;
        this.max_skipped_exon = builder.max_skipped_exon;
        this.min_skipped_bases = builder.min_skipped_bases;
        this.max_skipped_bases = builder.max_skipped_bases;
    }

    public static Set<ExonSkip> findExonSkippingEvents(GTFAnnotation gtfAnnotation) {
        Set<ExonSkip> exonSkips = new HashSet<>();
        gtfAnnotation.getGenes().values().parallelStream().forEach(gene -> {
            for (Transcript transcript : gene.getTranscripts().values()) {
                for (Interval intron : transcript.getIntrons()) {
                    // Check if the intron is a candidate for exon skipping
                    Set<String> spliceVariantProteins = new HashSet<>();
                    Set<String> wildTypeStartProteins = new HashSet<>();
                    Set<String> wildTypeEndProteins = new HashSet<>();
                    Set<Interval> wildTypeIntrons = new HashSet<>();

                    // TODO check if the intro has already been processed in a set
                    findMatchingTranscriptsForIntron(gene, intron, spliceVariantProteins, wildTypeStartProteins, wildTypeIntrons, wildTypeEndProteins);

                    // (WT_start intersect WT_end) – SV Folie 22/26
                    Set<String> wildTypeProteins = new HashSet<>(wildTypeStartProteins);
                    wildTypeProteins.retainAll(wildTypeEndProteins);
                    wildTypeProteins.removeAll(spliceVariantProteins);

                    if (!wildTypeProteins.isEmpty()) {
                        ExonSkip exonSkip = createExonSkipEvent(gene, intron, wildTypeProteins, wildTypeIntrons, spliceVariantProteins);
                        exonSkips.add(exonSkip);
                    }
                }
            }
        });
        return exonSkips;
    }

    private static ExonSkip createExonSkipEvent(Gene gene, Interval intron, Set<String> WT, Set<Interval> wildTypeIntrons, Set<String> spliceVariantProteins) {
        SkippedExonsBases skippedExonsBases = calculateSkippedExonsAndBases(gene, intron, WT);
        return new Builder()
                .setId(gene.getId())
                .setSymbol(gene.getAttribute(GTFParser.GENE_NAME))
                .setChr(gene.getSeqname())
                .setStrand(gene.getStrand())
                .setNprots((int) gene.getTranscripts().values().stream().filter(t -> !t.getCds().isEmpty()).count())
                .setNtrans(gene.getTranscripts().size())
                .setSV(new Interval(intron))
                .setWT(new TreeSet<>(wildTypeIntrons))
                .setSV_prots(spliceVariantProteins)
                .setWT_prots(WT)
                .setMax_skipped_exon(skippedExonsBases.maxSkippedExons())
                .setMin_skipped_exon(skippedExonsBases.minSkippedExons())
                .setMax_skipped_bases(skippedExonsBases.maxSkippedBases())
                .setMin_skipped_bases(skippedExonsBases.minSkippedBases())
                .build();
    }

    private static SkippedExonsBases calculateSkippedExonsAndBases(Gene gene, Interval intron, Set<String> WT) {
        int minSkippedExons = Integer.MAX_VALUE;
        int maxSkippedExons = Integer.MIN_VALUE;
        int minSkippedBases = Integer.MAX_VALUE;
        int maxSkippedBases = Integer.MIN_VALUE;
        for (String protein_id : WT) {
            // Get the transcript
            String transcript_id = gene.convertProteinToTranscript(protein_id);
            Transcript t = gene.getTranscript(transcript_id);
            int skippedExons = 0;
            int skippedBases = 0;
            // Find the number of exons that lie within the intron
            for (Exon exon : t.getExons()) {
                if (intron.getStart() <= exon.getInterval().getStart() && intron.getEnd() >= exon.getInterval().getEnd()) {
                    skippedExons++;
                    skippedBases += exon.getInterval().getLength();
                }
            }
            minSkippedExons = Math.min(minSkippedExons, skippedExons);
            maxSkippedExons = Math.max(maxSkippedExons, skippedExons);
            minSkippedBases = Math.min(minSkippedBases, skippedBases);
            maxSkippedBases = Math.max(maxSkippedBases, skippedBases);

        }
        return new SkippedExonsBases(minSkippedExons, maxSkippedExons, minSkippedBases, maxSkippedBases);
    }

    private static void findMatchingTranscriptsForIntron(Gene gene, Interval intron, Set<String> spliceVariantProteins, Set<String> wildTypeStartProteins, Set<Interval> wildTypeIntrons, Set<String> wildTypeEndProteins) {
        for (Transcript transcriptToCheck : gene.getTranscripts().values()) {
            // only consider transcripts with cds
            if (transcriptToCheck.getCds().isEmpty()) {
                //System.out.println("No CDS found for transcript " + t.getId());
                continue;
            }

            for (Interval intronToCheck : transcriptToCheck.getIntrons()) {
                // Check if the intron is the same as the candidate intron (SV)

                // TODO: check if an intron is in the middle (so doesn't match start or end)
                // TODO: do we even need to do the set operations
                if (intronToCheck.equals(intron)) {
                    spliceVariantProteins.add(transcriptToCheck.getCds().getFirst().getAttribute(GTFParser.PROTEIN_ID));
                } else if (intronToCheck.getStart() == intron.getStart()) {
                    wildTypeStartProteins.add(transcriptToCheck.getCds().getFirst().getAttribute(GTFParser.PROTEIN_ID));
                    wildTypeIntrons.add(intronToCheck);
                } else if (intronToCheck.getEnd() == intron.getEnd()) {
                    wildTypeEndProteins.add(transcriptToCheck.getCds().getFirst().getAttribute(GTFParser.PROTEIN_ID));
                    wildTypeIntrons.add(intronToCheck);
                    // if it is the end intron, add it to the wildTypeIntrons
                } else if (intronToCheck.getStart() > intron.getStart() && intronToCheck.getEnd() < intron.getEnd()) {
                    wildTypeIntrons.add(intronToCheck);
                }
            }
        }
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ExonSkip exonSkip = (ExonSkip) o;
        return nprots == exonSkip.nprots && ntrans == exonSkip.ntrans && min_skipped_exon == exonSkip.min_skipped_exon && max_skipped_exon == exonSkip.max_skipped_exon && min_skipped_bases == exonSkip.min_skipped_bases && max_skipped_bases == exonSkip.max_skipped_bases && Objects.equals(id, exonSkip.id) && Objects.equals(symbol, exonSkip.symbol) && Objects.equals(chr, exonSkip.chr) && strand == exonSkip.strand && Objects.equals(SV, exonSkip.SV) && Objects.equals(WT, exonSkip.WT) && Objects.equals(SV_prots, exonSkip.SV_prots) && Objects.equals(WT_prots, exonSkip.WT_prots);
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public int getMax_skipped_bases() {
        return max_skipped_bases;
    }

    public void setMax_skipped_bases(int max_skipped_bases) {
        this.max_skipped_bases = max_skipped_bases;
    }

    public int getMax_skipped_exon() {
        return max_skipped_exon;
    }

    public void setMax_skipped_exon(int max_skipped_exon) {
        this.max_skipped_exon = max_skipped_exon;
    }

    public int getMin_skipped_bases() {
        return min_skipped_bases;
    }

    public void setMin_skipped_bases(int min_skipped_bases) {
        this.min_skipped_bases = min_skipped_bases;
    }

    public int getMin_skipped_exon() {
        return min_skipped_exon;
    }

    public void setMin_skipped_exon(int min_skipped_exon) {
        this.min_skipped_exon = min_skipped_exon;
    }

    public int getNprots() {
        return nprots;
    }

    public void setNprots(int nprots) {
        this.nprots = nprots;
    }

    public int getNtrans() {
        return ntrans;
    }

    public void setNtrans(int ntrans) {
        this.ntrans = ntrans;
    }

    public StrandDirection getStrand() {
        return strand;
    }

    public void setStrand(StrandDirection strand) {
        this.strand = strand;
    }

    public Interval getSV() {
        return SV;
    }

    public void setSV(Interval SV) {
        this.SV = SV;
    }

    public Set<String> getSV_prots() {
        return SV_prots;
    }

    public void setSV_prots(Set<String> SV_prots) {
        this.SV_prots = SV_prots;
    }

    public String getSymbol() {
        return symbol;
    }

    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    public TreeSet<Interval> getWT() {
        return WT;
    }

    public void setWT(TreeSet<Interval> WT) {
        this.WT = WT;
    }

    public Set<String> getWT_prots() {
        return WT_prots;
    }

    public void setWT_prots(Set<String> WT_prots) {
        this.WT_prots = WT_prots;
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, symbol, chr, strand, nprots, ntrans, SV, WT, SV_prots, WT_prots, min_skipped_exon, max_skipped_exon, min_skipped_bases, max_skipped_bases);
    }

    private record SkippedExonsBases(int minSkippedExons, int maxSkippedExons, int minSkippedBases,
                                     int maxSkippedBases) {
    }

    public static class Builder {
        private String id;
        private String symbol;
        private String chr;
        private StrandDirection strand;
        private int nprots;
        private int ntrans;
        private Interval SV;
        private TreeSet<Interval> WT;
        private Set<String> SV_prots = new HashSet<>();
        private Set<String> WT_prots = new HashSet<>();
        private int min_skipped_exon;
        private int max_skipped_exon;
        private int min_skipped_bases;
        private int max_skipped_bases;

        public Builder setId(String id) {
            this.id = id;
            return this;
        }

        public Builder setSymbol(String symbol) {
            this.symbol = symbol;
            return this;
        }

        public Builder setChr(String chr) {
            this.chr = chr;
            return this;
        }

        public Builder setStrand(StrandDirection strand) {
            this.strand = strand;
            return this;
        }

        public Builder setNprots(int nprots) {
            this.nprots = nprots;
            return this;
        }

        public Builder setNtrans(int ntrans) {
            this.ntrans = ntrans;
            return this;
        }

        public Builder setSV(Interval SV) {
            this.SV = SV;
            return this;
        }

        public Builder setWT(TreeSet<Interval> WT) {
            this.WT = WT;
            return this;
        }

        public Builder setSV_prots(Set<String> SV_prots) {
            this.SV_prots = SV_prots;
            return this;
        }

        public Builder setWT_prots(Set<String> WT_prots) {
            this.WT_prots = WT_prots;
            return this;
        }

        public Builder setMin_skipped_exon(int min_skipped_exon) {
            this.min_skipped_exon = min_skipped_exon;
            return this;
        }

        public Builder setMax_skipped_exon(int max_skipped_exon) {
            this.max_skipped_exon = max_skipped_exon;
            return this;
        }

        public Builder setMin_skipped_bases(int min_skipped_bases) {
            this.min_skipped_bases = min_skipped_bases;
            return this;
        }

        public Builder setMax_skipped_bases(int max_skipped_bases) {
            this.max_skipped_bases = max_skipped_bases;
            return this;
        }

        public ExonSkip build() {
            return new ExonSkip(this);
        }
    }
}
