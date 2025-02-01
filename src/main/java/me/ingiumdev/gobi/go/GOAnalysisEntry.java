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

    public void setKs_fdr(double ks_fdr) {
        this.ks_fdr = ks_fdr;
    }

    public void setFej_fdr(double fej_fdr) {
        this.fej_fdr = fej_fdr;
    }

    public void setHg_fdr(double hg_fdr) {
        this.hg_fdr = hg_fdr;
    }

    public double getHg_pvalue() {
        return hg_pvalue;
    }

    public double getFej_pvalue() {
        return fej_pvalue;
    }

    public double getKs_pvalue() {
        return ks_pvalue;
    }

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
        return "GO:" + rawID + "\t" + name + "\t" + size + "\t" + isSOT
                + "\t" + numOverlap + "\t" + roundToFiveDecimalPlaces(hg_pvalue) + "\t" + roundToFiveDecimalPlaces(hg_fdr) + "\t" + roundToFiveDecimalPlaces(fej_pvalue) + "\t" + roundToFiveDecimalPlaces(fej_fdr) + "\t" + roundToFiveDecimalPlaces(ks_stat) + "\t" + roundToFiveDecimalPlaces(ks_pvalue) + "\t" + ks_fdr
                + "\t" + shortest_path_to_a_true
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
        // If no true entries are provided or if the analyzed term is already true, nothing to do.
        if (soTterms == null || soTterms.isEmpty() || soTterms.contains(goTerm.getId())) {
            this.shortest_path_to_a_true = "";
            return;
        }

        // 1. Upward search: record all nodes reachable by going upward (using only parent pointers)
        //    from the analyzed term. For each node, store the number of upward moves and a pointer
        //    to reconstruct the upward path.
        Map<GOTerm, Integer> upwardDist = new HashMap<>();
        Map<GOTerm, GOTerm> upwardPred = new HashMap<>();
        Queue<GOTerm> upwardQueue = new LinkedList<>();

        upwardQueue.add(goTerm);
        upwardDist.put(goTerm, 0);
        upwardPred.put(goTerm, null);

        while (!upwardQueue.isEmpty()) {
            GOTerm current = upwardQueue.poll();
            for (GOTerm parent : current.getParents()) {
                if (!upwardDist.containsKey(parent)) {
                    upwardDist.put(parent, upwardDist.get(current) + 1);
                    upwardPred.put(parent, current);
                    upwardQueue.add(parent);
                }
            }
        }

        // 2. For each candidate upward node (ancestor), attempt a downward search
        //    (using only child pointers) to reach a true node.
        //    We will record the total valid path length (upward moves + downward moves)
        //    and keep track of the best candidate (the turning point or LCA).
        int bestTotalDist = Integer.MAX_VALUE;
        GOTerm bestCandidate = null;
        GOTerm bestTrueNode = null;
        // For the best candidate, store its downward predecessor mapping for path reconstruction.
        Map<GOTerm, GOTerm> bestDownwardPred = null;

        for (GOTerm candidate : upwardDist.keySet()) {
            // Downward BFS from candidate (using only children pointers)
            Map<GOTerm, Integer> downwardDist = new HashMap<>();
            Map<GOTerm, GOTerm> downwardPred = new HashMap<>();
            Queue<GOTerm> downwardQueue = new LinkedList<>();

            downwardQueue.add(candidate);
            downwardDist.put(candidate, 0);
            downwardPred.put(candidate, null);

            GOTerm foundTrue = null;
            while (!downwardQueue.isEmpty()) {
                GOTerm current = downwardQueue.poll();
                if (soTterms.contains(current.getId())) {
                    foundTrue = current;
                    break;
                }
                for (GOTerm child : current.getChildren()) {
                    if (!downwardDist.containsKey(child)) {
                        downwardDist.put(child, downwardDist.get(current) + 1);
                        downwardPred.put(child, current);
                        downwardQueue.add(child);
                    }
                }
            }

            if (foundTrue != null) {
                int totalDist = upwardDist.get(candidate) + downwardDist.get(foundTrue);
                if (totalDist < bestTotalDist) {
                    bestTotalDist = totalDist;
                    bestCandidate = candidate;
                    bestTrueNode = foundTrue;
                    bestDownwardPred = downwardPred;
                }
            }
        }

        // If no candidate yields a valid downward path, then leave the path empty.
        if (bestCandidate == null || bestTrueNode == null) {
            this.shortest_path_to_a_true = "";
            return;
        }

        // 3. Reconstruct the upward path from the analyzed term to the best candidate (LCA).
        List<GOTerm> upwardPath = new ArrayList<>();
        GOTerm curr = bestCandidate;
        while (curr != null) {
            upwardPath.add(curr);
            curr = upwardPred.get(curr);
        }
        Collections.reverse(upwardPath); // Now, upwardPath[0] is goTerm, last element is the LCA.

        // 4. Reconstruct the downward path from the candidate (LCA) to the true term,
        //    using the stored downward predecessor mapping.
        List<GOTerm> downwardPath = new ArrayList<>();
        curr = bestTrueNode;
        while (curr != null && !curr.equals(bestCandidate)) {
            downwardPath.add(curr);
            curr = bestDownwardPred.get(curr);
        }
        Collections.reverse(downwardPath);
        // Avoid duplicating the candidate if it appears at the start of the downward path.
        if (!downwardPath.isEmpty() && downwardPath.get(0).equals(bestCandidate)) {
            downwardPath.remove(0);
        }

        // 5. Assemble the final path: join the upward part and the downward part,
        //    marking the turning point (LCA) with an appended " *".
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (int i = 0; i < upwardPath.size(); i++) {
            if (!first) {
                sb.append("|");
            }
            first = false;
            sb.append(upwardPath.get(i).getName());
            if (i == upwardPath.size() - 1) {  // this is the candidate (LCA)
                sb.append(" *");
            }
        }
        for (GOTerm term : downwardPath) {
            sb.append("|").append(term.getName());
        }

        this.shortest_path_to_a_true = sb.toString();
    }
}
