package gtf;

import java.util.List;
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
}
