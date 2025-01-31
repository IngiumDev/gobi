package me.ingiumdev.gobi.go;

import java.util.Map;

public class GOAnalysisEntry {

    private final int id;
    private final String rawID;
    private final String name;
    private int size;
    private final boolean isSOT;
    private int numOverlap;
    private double hg_pvalue;
    private double hg_fdr;
    private double fej_pvalue;
    private double fej_fdr;
    private double ks_stat;
    private double ks_pvalue;
    private double ks_fdr;
    private String shortest_path_to_a_true;

    public GOAnalysisEntry(GOTerm term) {
        this.id = term.getId();
        this.rawID = term.getFullID();
        this.name = term.getName();
        this.isSOT = term.isSOT();
    }

    public void calculateSize(DAG graph, Map<String, GeneSetEnrichmentAnalysis.DifferentialExpressionRecord> differentialExpressionInput) {
        /* size: number of measured associated genes to the GO categories (i.e. the
number of gene ids both occurring in the file given by the enrich option and
associated to the GO entry by the provided mapping (see option mapping).*/
        this.size = 0;
        for (String gene : differentialExpressionInput.keySet()) {
            if (graph.getTerm(id).getAssociatedGenes().contains(gene)) {
                this.size++;
                if (differentialExpressionInput.get(gene).signif()) {
                    this.numOverlap++;
                }
            }
        }
    }


    @Override
    public String toString() {
        return "GO:" + rawID + "\t" + name + "\t" + size + "\t" + isSOT + "\t" + numOverlap //+ "\t" + hg_pvalue + "\t" + hg_fdr + "\t" + fej_pvalue + "\t" + fej_fdr + "\t" + ks_stat + "\t" + ks_pvalue + "\t" + ks_fdr + "\t" + shortest_path_to_a_true
                ;
    }

    public int getId() {
        return id;
    }
}
