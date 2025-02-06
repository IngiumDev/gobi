package me.ingiumdev.gobi.go;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import java.util.*;
import java.util.stream.Collectors;

public class GOAnalysisEntry {

    private final int id;
    private final String rawID;
    private final String name;
    private final boolean isSOT;
    private int size;
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

    public static double roundToFiveDecimalPlaces(double value) {
        if (value == 0.0) {
            return 0.0;
        }
        double sign = Math.signum(value);
        double absValue = Math.abs(value);
        double exponent = Math.floor(Math.log10(absValue));
        double scale = Math.pow(10, 5 - exponent);
        double rounded = Math.round(absValue * scale) / scale;
        return sign * rounded;
    }

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

    public int getId() {
        return id;
    }

    @Override
    public String toString() {
        return "GO:" + rawID + "\t" + name + "\t" + size + "\t" + isSOT + "\t" + numOverlap + "\t" + roundToFiveDecimalPlaces(hg_pvalue) + "\t" + roundToFiveDecimalPlaces(hg_fdr) + "\t" + roundToFiveDecimalPlaces(fej_pvalue) + "\t" + roundToFiveDecimalPlaces(fej_fdr) + "\t" + roundToFiveDecimalPlaces(ks_stat) + "\t" + roundToFiveDecimalPlaces(ks_pvalue) + "\t" + ks_fdr + "\t" + shortest_path_to_a_true;
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

        List<Double> measuredGeneFC = sizeGenes.stream().map(gene -> differentialExpressionInput.get(gene).fc()).collect(Collectors.toCollection(() -> new ArrayList<>(sizeGenes.size())));
        List<Double> backgroundGeneFC = graph.getRoot().getAssociatedGenes().stream().filter(t -> !sizeGenes.contains(t) && differentialExpressionInput.containsKey(t)).map(t -> differentialExpressionInput.get(t).fc()).collect(Collectors.toCollection(() -> new ArrayList<>(graph.getRoot().getAssociatedGenes().size())));

        double[] measuredArray = new double[measuredGeneFC.size()];
        for (int i = 0; i < measuredGeneFC.size(); i++) {
            measuredArray[i] = measuredGeneFC.get(i);
        }

        double[] backgroundArray = new double[backgroundGeneFC.size()];
        for (int i = 0; i < backgroundGeneFC.size(); i++) {
            backgroundArray[i] = backgroundGeneFC.get(i);
        }

        ks_stat = ks.kolmogorovSmirnovStatistic(measuredArray, backgroundArray);
        ks_pvalue = ks.kolmogorovSmirnovTest(measuredArray, backgroundArray);
    }


    public void calculateShortestPathToTrue(Set<Integer> SOTTerms, GOTerm goTerm) {
        // If no true entries are provided or if the analyzed term is already true, nothing to do.
        if (SOTTerms == null || SOTTerms.isEmpty() || SOTTerms.contains(goTerm.getId())) {
            this.shortest_path_to_a_true = "";
            return;
        }

        // 1. Upward search: record all nodes reachable by going upward from the analyzed term.
        // For each node, store the number of upward moves and a pointer
        // to reconstruct the upward path.
        // upwardDist: distance from the analyzed term to the current node.
        // upwardPred: predecessor node in the upward path.
        Map<GOTerm, Integer> upwardDist = new HashMap<>();
        Map<GOTerm, GOTerm> upwardPred = new HashMap<>();
        // Use ArrayDeque for faster queue operations.
        // --> special kind of array that grows
        // and allows users to add or remove an element from both sides of the queue.
        // Queue stores the nodes to be processed.
        Queue<GOTerm> upwardQueue = new ArrayDeque<>();

        upwardQueue.add(goTerm);
        upwardDist.put(goTerm, 0);
        upwardPred.put(goTerm, null);

        while (!upwardQueue.isEmpty()) {
            GOTerm current = upwardQueue.poll();
            int currentDist = upwardDist.get(current);
            for (GOTerm parent : current.getParents()) {
                if (!upwardDist.containsKey(parent)) {
                    upwardDist.put(parent, currentDist + 1);
                    upwardPred.put(parent, current);
                    upwardQueue.add(parent);
                }
            }
        }

        // 2. Process candidates in increasing order of upward distance to allow for pruning.
        List<GOTerm> candidates = new ArrayList<>(upwardDist.keySet());
        candidates.sort(Comparator.comparingInt(upwardDist::get));

        int bestTotalDist = Integer.MAX_VALUE;
        GOTerm bestCandidate = null;
        GOTerm bestTrueNode = null;
        Map<GOTerm, GOTerm> bestDownwardPred = null;

        // For each candidate (potential turning point / LCA): check if there is a downward path to a true node.
        for (GOTerm candidate : candidates) {
            int candidateUpwardDist = upwardDist.get(candidate);
            // Prune candidates that already have an upward cost not lower than the current best.
            if (candidateUpwardDist >= bestTotalDist) {
                break;
            }

            // Downward BFS from candidate (using only child pointers).
            Map<GOTerm, Integer> downwardDist = new HashMap<>();
            Map<GOTerm, GOTerm> downwardPred = new HashMap<>();
            Queue<GOTerm> downwardQueue = new ArrayDeque<>();

            downwardQueue.add(candidate);
            downwardDist.put(candidate, 0);
            downwardPred.put(candidate, null);

            GOTerm foundTrue = null;
            while (!downwardQueue.isEmpty()) {
                GOTerm current = downwardQueue.poll();
                int currentDownwardDist = downwardDist.get(current);
                // Prune: if the sum already meets or exceeds bestTotalDist, skip further processing from this node.
                if (candidateUpwardDist + currentDownwardDist >= bestTotalDist) {
                    continue;
                }
                if (SOTTerms.contains(current.getId())) {
                    foundTrue = current;
                    break;
                }
                for (GOTerm child : current.getChildren()) {
                    if (!downwardDist.containsKey(child)) {
                        int newDownwardDist = currentDownwardDist + 1;
                        // Prune early: if candidateUpwardDist + newDownwardDist is already worse, skip this child.
                        if (candidateUpwardDist + newDownwardDist >= bestTotalDist) {
                            continue;
                        }
                        downwardDist.put(child, newDownwardDist);
                        downwardPred.put(child, current);
                        downwardQueue.add(child);
                    }
                }
            }

            if (foundTrue != null) {
                int totalDist = candidateUpwardDist + downwardDist.get(foundTrue);
                if (totalDist < bestTotalDist) {
                    bestTotalDist = totalDist;
                    bestCandidate = candidate;
                    bestTrueNode = foundTrue;
                    bestDownwardPred = downwardPred;
                }
            }
        }

        // If no candidate yields a valid downward path, then leave the path empty.
        if (bestCandidate == null) {
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

        // 4. Reconstruct the downward path from the candidate (LCA) to the true node,
        //    using the stored downward predecessor mapping.
        List<GOTerm> downwardPath = new ArrayList<>();
        curr = bestTrueNode;
        while (curr != null && !curr.equals(bestCandidate)) {
            downwardPath.add(curr);
            curr = bestDownwardPred.get(curr);
        }
        Collections.reverse(downwardPath);
        // Avoid duplicating the candidate if it appears at the start of the downward path.
        if (!downwardPath.isEmpty() && downwardPath.getFirst().equals(bestCandidate)) {
            downwardPath.removeFirst();
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
