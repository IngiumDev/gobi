package gtf;

import java.util.List;
import java.util.Objects;
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
    private List<String> SV_prots;
    // ids of the WT CDS, separated by | (protein_id)
    private List<String> WT_prots;
    // the minimal number of skipped exons in any WT/SV pair
    private int min_skipped_exon;
    // the maximum number of skipped exons in any WT/SV pair
    private int max_skipped_exon;
    // the minimal number of skipped bases (joint length of skipped exons) in any WT/SV pair
    private int min_skipped_bases;
    // the maximum number of skipped bases (joint length of skipped exons) in any WT/SV pair
    private int max_skipped_bases;

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

    public List<String> getSV_prots() {
        return SV_prots;
    }

    public void setSV_prots(List<String> SV_prots) {
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

    public List<String> getWT_prots() {
        return WT_prots;
    }

    public void setWT_prots(List<String> WT_prots) {
        this.WT_prots = WT_prots;
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, symbol, chr, strand, nprots, ntrans, SV, WT, SV_prots, WT_prots, min_skipped_exon, max_skipped_exon, min_skipped_bases, max_skipped_bases);
    }
}
