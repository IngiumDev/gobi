package gtf;

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

    public static Set<ExonSkip> findExonSkippingEvents(GTFAnnotation GTFAnnotation) {
        Set<ExonSkip> exonSkips = new HashSet<>();
        GTFAnnotation.getGenes().values().parallelStream().forEach(gene -> {
            for (Transcript transcript : gene.getTranscripts().values()) {
                for (Interval intron : transcript.getIntrons()) {
                    // Check if the intron is a candidate for exon skipping
                    Set<String> SV = new HashSet<>();
                    Set<String> WT_start = new HashSet<>();
                    Set<String> WT_end = new HashSet<>();
                    Set<Interval> WT_introns = new HashSet<>();

                    int minSkippedExons = Integer.MAX_VALUE;
                    int maxSkippedExons = Integer.MIN_VALUE;
                    int minSkippedBases = Integer.MAX_VALUE;
                    int maxSkippedBases = Integer.MIN_VALUE;

                    for (Transcript t : gene.getTranscripts().values()) {
                        // only consider transcripts with cds
                        if (t.getCds().isEmpty()) {
                            //System.out.println("No CDS found for transcript " + t.getId());
                            continue;
                        }

                        for (Interval i : t.getIntrons()) {
                            // Check if the intron is the same as the candidate intron (SV)
                            String protein_id;
                            if (i.equals(intron)) {
                                protein_id = t.getCds().getFirst().getAttribute("protein_id");
                                SV.add(protein_id);
                            } else if (i.getStart() == intron.getStart()) {
                                protein_id = t.getCds().getFirst().getAttribute("protein_id");
                                WT_start.add(protein_id);
                                WT_introns.add(i);
                            } else if (i.getEnd() == intron.getEnd()) {
                                protein_id = t.getCds().getFirst().getAttribute("protein_id");
                                WT_end.add(protein_id);
                                WT_introns.add(i);
                            }
                        }
                    }
                    Set<String> WT = new HashSet<>(WT_start);
                    WT.retainAll(WT_end);
                    WT.removeAll(SV);
                    if (!WT.isEmpty()) {
                        // We have an exon skipping event
                        // TODO add exon skip to the list
                        ExonSkip exonSkip = new ExonSkip();
                        exonSkip.setId(gene.getId());
                        // If gene is directly annoted
                        exonSkip.setSymbol(gene.getAttribute("gene_name"));
                        exonSkip.setChr(gene.getSeqname());
                        exonSkip.setStrand(gene.getStrand());

                        // nprots is the number of transcripts where the cds is not empty
                        exonSkip.setNprots((int) gene.getTranscripts().values().stream().filter(t -> !t.getCds().isEmpty()).count());
                        exonSkip.setNtrans(gene.getTranscripts().size());
                        exonSkip.setSV(new Interval(intron));
                        exonSkip.setWT(new TreeSet<>(WT_introns));
                        exonSkip.setSV_prots(SV);
                        exonSkip.setWT_prots(WT);
                        // Calculate min and max skipped exons and bases

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

                        exonSkip.setMax_skipped_exon(maxSkippedExons);
                        exonSkip.setMin_skipped_exon(minSkippedExons);
                        exonSkip.setMax_skipped_bases(maxSkippedBases);
                        exonSkip.setMin_skipped_bases(minSkippedBases);

                        exonSkips.add(exonSkip);
                    }
                }
            }
        });
        return exonSkips;
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
}
