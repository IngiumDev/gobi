package gtf;

import gtf.structs.CodingSequence;
import gtf.structs.Gene;
import gtf.structs.Interval;
import gtf.structs.Transcript;
import gtf.types.StrandDirection;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

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
        long start = System.currentTimeMillis();
        gtfAnnotation.getGenes().values().parallelStream().forEach(gene -> {
            for (Interval intron : gene.getIntrons()) {
                Set<String> spliceVariantTranscripts = new HashSet<>();
                Set<String> wildTypeTranscripts = new HashSet<>();
                Set<String> wildTypeEndTranscripts = new HashSet<>();
                Set<Interval> wildTypeIntrons = new TreeSet<>();
                checkIntronCandidate(gene, intron, spliceVariantTranscripts, wildTypeTranscripts, wildTypeIntrons);

                if (!spliceVariantTranscripts.isEmpty() && !wildTypeTranscripts.isEmpty()) {
                    ExonSkip exonSkip = createExonSkipEvent(gene, intron, spliceVariantTranscripts, wildTypeIntrons, wildTypeTranscripts);
                    exonSkips.add(exonSkip);
                }
            }
        });
        System.out.println("Time taken: " + (System.currentTimeMillis() - start) + "ms");
        return exonSkips;
    }

    private static void checkIntronCandidate(Gene gene, Interval intronCandidate, Set<String> spliceVariantTranscripts, Set<String> wildTypeTranscripts, Set<Interval> wildTypeIntrons) {
        for (Transcript transcriptToCheck : gene.getTranscripts().values()) {
            Set<Interval> intronsToAdd = new TreeSet<>();

            int intronCandidateStart = intronCandidate.getStart();
            int intronCandidateEnd = intronCandidate.getEnd();
            boolean isStartIntron = false;
            boolean isEndIntron = false;
            CodingSequence previousCds = null;




            for (CodingSequence cds : transcriptToCheck.getCds()) {
                if (previousCds != null) {
                    if (previousCds.getInterval().getEnd() + 1 >= intronCandidateStart && cds.getInterval().getStart() - 1 <= intronCandidateEnd ) {
                        if (previousCds.getInterval().getEnd() + 1 == intronCandidateStart && cds.getInterval().getStart() -1 == intronCandidateEnd) {
                            spliceVariantTranscripts.add(transcriptToCheck.getTranscriptID());
                            break;
                        } else if (previousCds.getInterval().getEnd() + 1 == intronCandidateStart) {
                            isStartIntron = true;
                            intronsToAdd.add(new Interval(intronCandidateStart, cds.getInterval().getStart() - 1));
                        } else if (cds.getInterval().getStart() - 1 == intronCandidateEnd) {
                            isEndIntron = true;
                            intronsToAdd.add(new Interval(previousCds.getInterval().getEnd() +1, intronCandidateEnd));
                        } else {
                            intronsToAdd.add(new Interval(previousCds.getInterval().getEnd() + 1, cds.getInterval().getStart() - 1));
                        }
                    }
                }
                previousCds = cds;
            }
            if (isStartIntron && isEndIntron) {
                wildTypeTranscripts.add(transcriptToCheck.getTranscriptID());
                wildTypeIntrons.addAll(intronsToAdd);
            }
        }

    }

    private static ExonSkip createExonSkipEvent(Gene gene, Interval intron, Set<String> wildTypeTranscripts, Set<Interval> wildTypeIntrons, Set<String> spliceVariantTranscripts) {
        SkippedExonsBases skippedExonsBases = calculateSkippedExonsAndBases(gene, intron, spliceVariantTranscripts);
        Set<String> spliceVariantProteins = convertTranscriptToProteinID(gene, spliceVariantTranscripts);
        Set<String> wildTypeProteins = convertTranscriptToProteinID(gene, wildTypeTranscripts);
        return new Builder()
                .setId(gene.getGeneID())
                .setSymbol(gene.getGeneName())
                .setChr(gene.getSeqname())
                .setStrand(gene.getStrand())
                .setNprots((int) gene.getTranscripts().values().stream().map(Transcript::getCds).filter(cds -> cds != null && !cds.isEmpty()).count())
                .setNtrans(gene.getTranscripts().size())
                .setSV(new Interval(intron))
                .setWT(new TreeSet<>(wildTypeIntrons))
                .setSV_prots(spliceVariantProteins)
                .setWT_prots(wildTypeProteins)
                .setMax_skipped_exon(skippedExonsBases.maxSkippedExons())
                .setMin_skipped_exon(skippedExonsBases.minSkippedExons())
                .setMax_skipped_bases(skippedExonsBases.maxSkippedBases())
                .setMin_skipped_bases(skippedExonsBases.minSkippedBases())
                .build();
    }

    private static Set<String> convertTranscriptToProteinID(Gene gene, Set<String> transcriptIDs) {
        return transcriptIDs.stream()
                .map(gene::getTranscript)
                .map(transcript -> {
                    String proteinID = transcript.getCds().getFirst().getProteinID();
                    if (proteinID != null) {
                        return proteinID;
                    }
                    String ccdsID = transcript.getCds().getFirst().getCcdsID();
                     if (ccdsID != null) {
                         return ccdsID;
                     }                      return "NONE";
                }).collect(Collectors.toSet());
    }

    private static SkippedExonsBases calculateSkippedExonsAndBases(Gene gene, Interval intron, Set<String> wildTypeTranscripts) {
        int minSkippedExons = Integer.MAX_VALUE;
        int maxSkippedExons = Integer.MIN_VALUE;
        int minSkippedBases = Integer.MAX_VALUE;
        int maxSkippedBases = Integer.MIN_VALUE;
        for (String transcript_id : wildTypeTranscripts) {
            // Get the transcript
            Transcript t = gene.getTranscript(transcript_id);
            int skippedExons = 0;
            int skippedBases = 0;
            // Find the number of exons that lie within the intron
            for (CodingSequence cds : t.getCds()) {

                if (intron.getStart() <= cds.getInterval().getStart() && intron.getEnd() >= cds.getInterval().getEnd()) {
                    skippedExons++;
                    skippedBases += cds.getInterval().getLength();
                }
            }
            minSkippedExons = Math.min(minSkippedExons, skippedExons);
            maxSkippedExons = Math.max(maxSkippedExons, skippedExons);
            minSkippedBases = Math.min(minSkippedBases, skippedBases);
            maxSkippedBases = Math.max(maxSkippedBases, skippedBases);

        }
        return new SkippedExonsBases(minSkippedExons, maxSkippedExons, minSkippedBases, maxSkippedBases);
    }

   /* private static void checkIntronCandidate(Gene gene, Interval intronCandidate, Set<String> spliceVariantProteins, Set<String> wildTypeStartProteins, Set<Interval> wildTypeIntrons, Set<String> wildTypeEndProteins) {
        for (Transcript transcriptToCheck : gene.getTranscripts().values()) {
            // only consider transcripts with cds
            if (transcriptToCheck.getCds().isEmpty()) {
                //System.out.println("No CDS found for transcript " + t.getId());
                continue;
            }

            for (Interval intronToCheck : transcriptToCheck.getIntrons()) {
                // Check if the intron is the same as the candidate intron (SV)

                if (intronToCheck.equals(intronCandidate)) {
                    spliceVariantProteins.add(transcriptToCheck.getCds().getFirst().getAttribute(GTFParser.PROTEIN_ID));
                } else if (intronToCheck.getStart() == intronCandidate.getStart()) {
                    wildTypeStartProteins.add(transcriptToCheck.getCds().getFirst().getAttribute(GTFParser.PROTEIN_ID));
                    wildTypeIntrons.add(intronToCheck);
                } else if (intronToCheck.getEnd() == intronCandidate.getEnd()) {
                    wildTypeEndProteins.add(transcriptToCheck.getCds().getFirst().getAttribute(GTFParser.PROTEIN_ID));
                    wildTypeIntrons.add(intronToCheck);
                    // if it is the end intron, add it to the wildTypeIntrons
                } else if (intronToCheck.getStart() > intronCandidate.getStart() && intronToCheck.getEnd() < intronCandidate.getEnd()) {
                    wildTypeIntrons.add(intronToCheck);
                }
            }
        }
    }*/


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
