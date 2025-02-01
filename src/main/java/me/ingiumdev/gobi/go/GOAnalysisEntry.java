package me.ingiumdev.gobi.go;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import java.util.*;

public class GOAnalysisEntry {

    private final int id;
    private final String rawID;
    private final String name;
    private int size;
    private final boolean isSOT;
    private int numOverlap;
    private Set<String> sizeGenes;
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
        sizeGenes = new HashSet<>();
        for (String gene : differentialExpressionInput.keySet()) {
            if (graph.getTerm(id).getAssociatedGenes().contains(gene)) {
                this.size++;
                sizeGenes.add(gene);
                if (differentialExpressionInput.get(gene).signif()) {
                    this.numOverlap++;
                }
            }
        }
    }


    public static double roundToFiveDecimalPlaces(double value) {
        if (value == 0.0) {
            return 0.0;
        }

        // Determine the sign of the value
        double sign = Math.signum(value);
        double absValue = Math.abs(value);

        // Calculate the base-10 exponent
        double exponent = Math.floor(Math.log10(absValue));

        // Scale the value to have five decimal places in the mantissa
        double scale = Math.pow(10, 5 - exponent);

        // Perform rounding
        double rounded = Math.round(absValue * scale) / scale;

        // Restore the original sign
        return sign * rounded;
    }

    public int getId() {
        return id;
    }

    @Override
    public String toString() {
        return "GO:" + rawID + "\t" + name + "\t" + size + "\t" + isSOT + "\t" + numOverlap + "\t" + roundToFiveDecimalPlaces(hg_pvalue) + "\t" + roundToFiveDecimalPlaces(hg_fdr) + "\t" + roundToFiveDecimalPlaces(fej_pvalue) + "\t" + fej_fdr + "\t" + roundToFiveDecimalPlaces(ks_stat) + "\t" + roundToFiveDecimalPlaces(ks_pvalue)// + "\t" + ks_fdr + "\t" + shortest_path_to_a_true
                ;
    }

    public void calculatehg_pvalue(int numDiffExpressedRootGenes, int numGenesInDifferentialExpression, GOTerm goTerm) {
        hg_pvalue = calculateHypergeometricPValue(numDiffExpressedRootGenes, numGenesInDifferentialExpression, size, numOverlap);
    }

    public void calculatefej_pvalue(int numDiffExpressedRootGenes, int numGenesInDifferentialExpression, GOTerm goTerm) {
        fej_pvalue = calculateHypergeometricPValue(numDiffExpressedRootGenes - 1, numGenesInDifferentialExpression - 1, size - 1, numOverlap - 1);
    }

    public double calculateHypergeometricPValue(int N, int K, int n, int k) {
        // N, K, n, k
        HypergeometricDistribution hypergeometricDistribution = new HypergeometricDistribution(N, K, n);
        return hypergeometricDistribution.upperCumulativeProbability(k);

    }

    public void calculateKS(DAG graph, Map<String, GeneSetEnrichmentAnalysis.DifferentialExpressionRecord> differentialExpressionInput) {
        KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
        List<Double> measuredGeneFC = new ArrayList<>();
        List<Double> backgroundGeneFC = new ArrayList<>();
        for (String t : graph.getRoot().getAssociatedGenes()) {
            if (!sizeGenes.contains(t) && differentialExpressionInput.containsKey(t)) {
                backgroundGeneFC.add(differentialExpressionInput.get(t).fc());
            }
        }
        for (var gene : sizeGenes) {
            measuredGeneFC.add(differentialExpressionInput.get(gene).fc());
        }
        ks_stat = ks.kolmogorovSmirnovStatistic(measuredGeneFC.stream().mapToDouble(Double::doubleValue).toArray(), backgroundGeneFC.stream().mapToDouble(Double::doubleValue).toArray());
        ks_pvalue = ks.kolmogorovSmirnovTest(measuredGeneFC.stream().mapToDouble(Double::doubleValue).toArray(), backgroundGeneFC.stream().mapToDouble(Double::doubleValue).toArray());
    }


    public void calculateShortestPathToTrue(Set<Integer> soTterms, DAG graph, GOTerm goTerm) {
        if (!soTterms.isEmpty() && !soTterms.contains(goTerm.getId())) {

        }
    }
}
