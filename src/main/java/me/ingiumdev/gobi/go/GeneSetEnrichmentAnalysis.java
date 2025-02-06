package me.ingiumdev.gobi.go;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class GeneSetEnrichmentAnalysis {
    private static final Logger log = LoggerFactory.getLogger(GeneSetEnrichmentAnalysis.class);
    private final ConcurrentMap<GOTerm, HashMap<GOTerm, Integer>> globalAncestorCache = new ConcurrentHashMap<>();
    int numDiffExpressedRootGenes;
    int numDiffExpressedSigRootGenes;
    private DAG graph;
    private Mapping mapping;
    private RootType rootType;
    private int minSize;
    private int maxSize;
    private Map<String, DifferentialExpressionRecord> differentialExpressionInput;
    private Set<Integer> SOTterms;
    private BufferedWriter writer;
    private String output;
    private String overlapOut;

    public GeneSetEnrichmentAnalysis(RootType rootType) {
    }


    private GeneSetEnrichmentAnalysis(Builder builder) {
        rootType = builder.rootType;
        minSize = builder.minSize;
        maxSize = builder.maxSize;
        output = builder.output;
        overlapOut = builder.overlapOut;
        try {
            writer = new BufferedWriter(new FileWriter(output));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static <T> void applyBenjaminiHochbergCorrection(List<T> entries, Function<T, Double> pValueGetter, BiConsumer<T, Double> fdrSetter) {
        int m = entries.size();

        // Sort the entries in ascending order based on the extracted p-value.
        List<T> sortedEntries = entries.stream().sorted(Comparator.comparingDouble(pValueGetter::apply)).toList();

        double[] adjustedFDR = new double[m];

        // Iterate backwards to compute the BH-adjusted p-values.
        for (int i = m - 1; i >= 0; i--) {
            double pValue = pValueGetter.apply(sortedEntries.get(i));
            int rank = i + 1; // Ranks are 1-based
            double bhValue = pValue * m / rank;

            if (i == m - 1) {
                adjustedFDR[i] = bhValue;
            } else {
                // take min
                adjustedFDR[i] = Math.min(bhValue, adjustedFDR[i + 1]);
            }
        }

        // Set the computed FDR values back to each entry.
        for (int i = 0; i < m; i++) {
            fdrSetter.accept(sortedEntries.get(i), adjustedFDR[i]);
        }
    }

    /**
     * Adjusts all the p-values for the provided GOAnalysisEntry list.
     *
     * @param analysisEntries The list of GOAnalysisEntry objects.
     */
    public static void correctPValues(List<GOAnalysisEntry> analysisEntries) {
        // Apply BH correction for the hypergeometric p-values.
        applyBenjaminiHochbergCorrection(analysisEntries, GOAnalysisEntry::getHg_pvalue, GOAnalysisEntry::setHg_fdr);

        // Apply BH correction for Fischer's exact test jackknife p-values.
        applyBenjaminiHochbergCorrection(analysisEntries, GOAnalysisEntry::getFej_pvalue, GOAnalysisEntry::setFej_fdr);

        // Apply BH correction for the KS test p-values.
        applyBenjaminiHochbergCorrection(analysisEntries, GOAnalysisEntry::getKs_pvalue, GOAnalysisEntry::setKs_fdr);
    }

    public void initMapping(String mappingPath, String mappingType) {
        long start = System.currentTimeMillis();
        if (mappingType.equals("ensembl")) {
            mapping = Mapping.createEnsemblMapping(mappingPath);
        } else {
            mapping = Mapping.createGOMapping(mappingPath);
        }
        log.info("Mapping loaded in {} ms", System.currentTimeMillis() - start);

    }

    public void initDAG(String oboPath) {
        graph = new DAG(oboPath, mapping, rootType);
    }

    public void initDifferentialExpression(String path) {
        long start = System.currentTimeMillis();
        Map<String, DifferentialExpressionRecord> diffExpRecords = new HashMap<>();
        SOTterms = new HashSet<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;

            while ((line = br.readLine()) != null) {
                int length = line.length();
                int i = 0;

                // Check if the line starts with #
                if (length > 0 && line.charAt(0) == '#') {
                    // Extract GO ID
                    i = 1; // Skip the #
                    while (i < length && line.charAt(i) == ' ') i++; // Skip leading spaces

                    int goStart = i;
                    while (i < length && line.charAt(i) != '\t' && line.charAt(i) != ' ') i++; // Find end of GO ID

                    if (goStart < i) {
                        String goID = line.substring(goStart, i);
                        if (goID.startsWith("GO:")) {
                            int goIntID = Integer.parseInt(goID.substring(3)); // Remove "GO:"
                            SOTterms.add(goIntID);
                            GOTerm goEntry = graph.getTerm(goIntID);
                            if (goEntry != null) {
                                goEntry.setSOT(); // Mark this GO term as Standard of Truth (SOT)
                            }
                        }
                    }
                    continue; // Skip to next line
                } else if (line.startsWith("id\tfc\tsignif")) {
                    continue;
                }

                // Manual parsing for normal data lines
                int idStart = i;
                while (i < length && line.charAt(i) != '\t') i++; // Find tab
                String geneID = line.substring(idStart, i);
                i++; // Skip tab

                // Parse FC (log2 fold change)
                int fcStart = i;
                while (i < length && line.charAt(i) != '\t') i++;
                double fc = Double.parseDouble(line.substring(fcStart, i));
                i++; // Skip tab

                // Parse significance (boolean)
                boolean signif = false;
                if (i < length) {
                    signif = line.charAt(i) == 't'; // Assume "true"/"false" values, so check first character
                }

                // Store the parsed record
                diffExpRecords.put(geneID, new DifferentialExpressionRecord(fc, signif));
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // Assign the parsed data to the class variable
        this.differentialExpressionInput = diffExpRecords;
        log.info("Loaded {} differential expression records in {} ms", diffExpRecords.size(), System.currentTimeMillis() - start);
    }

    public List<GOAnalysisEntry> performEnrichment() {
        long start = System.currentTimeMillis();
        calcNumDiffExpressedRootGenes(graph, differentialExpressionInput);
        List<GOAnalysisEntry> analysisEntries = graph.getEntries().values().parallelStream().filter(goTerm -> minSize <= goTerm.getAssociatedGenes().size() && goTerm.getAssociatedGenes().size() <= maxSize).map(this::analyzeGOEntry).collect(Collectors.toList());
        log.info("Calculated basic enrichment statistics in {} ms", System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        correctPValues(analysisEntries);
        log.info("Adjusted p-values in {} ms", System.currentTimeMillis() - start);
        return analysisEntries;
    }

    private void calcNumDiffExpressedRootGenes(DAG graph, Map<String, DifferentialExpressionRecord> differentialExpressionInput) {
        numDiffExpressedRootGenes = 0;
        numDiffExpressedSigRootGenes = 0;
        for (var differentialExpression : differentialExpressionInput.entrySet()) {
            if (graph.getRoot().getAssociatedGenes().contains(differentialExpression.getKey())) {
                numDiffExpressedRootGenes++;
                if (differentialExpression.getValue().signif()) {
                    numDiffExpressedSigRootGenes++;
                }
            }
        }

    }

    private GOAnalysisEntry analyzeGOEntry(GOTerm goTerm) {
        GOAnalysisEntry entry = new GOAnalysisEntry(goTerm);
        entry.calculateSize(graph, differentialExpressionInput);
        entry.calculatehg_pvalue(numDiffExpressedRootGenes, numDiffExpressedSigRootGenes, goTerm);
        entry.calculatefej_pvalue(numDiffExpressedRootGenes, numDiffExpressedSigRootGenes, goTerm);
        entry.calculateKS(graph, differentialExpressionInput);
        entry.calculateShortestPathToTrue(SOTterms, goTerm);
        return entry;
    }

    public void writeResults(List<GOAnalysisEntry> results) {
        long start = System.currentTimeMillis();
        try {
            writer.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true");
            for (GOAnalysisEntry entry : results) {
                writer.newLine();
                writer.write(entry.toString());
            }
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        log.info("Output written in {} ms", System.currentTimeMillis() - start);
    }

    public void performOverlapAnalysis() {
        long start = System.currentTimeMillis();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(overlapOut))) {
            writer.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent");
            // Retrieve and filter GOTerm objects based on the number of associated genes.
            List<GOTerm> filteredTerms = graph.getEntries().values().stream().filter(goTerm -> minSize <= goTerm.getAssociatedGenes().size() && goTerm.getAssociatedGenes().size() <= maxSize).toList();

            // Process unique pairs
            IntStream.range(0, filteredTerms.size()).parallel().forEach(i -> {
                for (int j = i + 1; j < filteredTerms.size(); j++) {
                    // Check if there is an overlap between the two GO terms.
                    int numSharedGenes = 0;
                    Set<String> associatedGenes = filteredTerms.get(i).getAssociatedGenes();
                    for (String gene : associatedGenes) {
                        if (filteredTerms.get(j).getAssociatedGenes().contains(gene)) {
                            numSharedGenes++;
                        }
                    }

                    if (numSharedGenes != 0) {
                        String result = processOverlapPair(filteredTerms.get(i), filteredTerms.get(j), numSharedGenes);
                        // Synchronized write to ensure thread safety.
                        synchronized (writer) {
                            try {
                                writer.write("\n" + result);
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        }
                    }
                }
            });
            System.out.println(filteredTerms.size());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        log.info("Overlap analysis written in {} ms", System.currentTimeMillis() - start);
    }

    /**
     * @param goTerm1
     * @param goTerm2
     * @param numSharedGenes
     * @return The principle is to restrict the search to the upward (parent) direction for each term, then combine the results.
     * For each term, you perform a directed search upward—tracking every ancestor and how many steps it takes to get there.
     * This gives you two maps: one for each term, with ancestors and their distances.
     * By intersecting these maps (i.e., finding common ancestors), you determine the shared nodes that lie in both terms’ ancestry chains.
     * For each common ancestor, adding the two distances yields a candidate path length (the sum of the upward moves from each term to that common ancestor).
     * The minimal sum over all common ancestors is taken as the shortest path length between the two terms.
     * This method is efficient because it only explores the typically small set of ancestors rather than the entire graph.
     * Comment by ChatGPT
     */
    private String processOverlapPair(GOTerm goTerm1, GOTerm goTerm2, int numSharedGenes) {
        int pathLength = computeShortestPath(goTerm1, goTerm2);
        boolean isRel = isRelative(goTerm1, goTerm2);
        double maxOvPercent = (double) numSharedGenes / Math.min(goTerm1.getAssociatedGenes().size(), goTerm2.getAssociatedGenes().size()) * 100;
        return "GO:" + goTerm1.getFullID() + "\t" + "GO:" + goTerm2.getFullID() + "\t" + isRel + "\t" + pathLength + "\t" + numSharedGenes + "\t" + maxOvPercent;
    }

    /**
     * @param term
     * @return returns a map from each visited ancestor to its distance (number of upward moves).
     * Computes and caches the set of ancestors for a given term (using only parent links)
     * and returns a map from each visited ancestor to its distance (number of upward moves)
     * from the start term. The cache is used only within one call to processOverlapPair.
     */
    private HashMap<GOTerm, Integer> computeAncestors(GOTerm term) {
        // computeIfAbsent ensures that for a given term only one thread computes its ancestors.
        return globalAncestorCache.computeIfAbsent(term, t -> {
            HashMap<GOTerm, Integer> distanceMap = new HashMap<>();
            Deque<GOTerm> queue = new ArrayDeque<>();
            distanceMap.put(t, 0);
            queue.add(t);
            while (!queue.isEmpty()) {
                GOTerm current = queue.poll();
                int currentDistance = distanceMap.get(current);
                for (GOTerm parent : current.getParents()) {
                    if (!distanceMap.containsKey(parent)) {
                        distanceMap.put(parent, currentDistance + 1);
                        queue.add(parent);
                    }
                }
            }
            return distanceMap;
        });
    }

    /**
     * @param term1
     * @param term2
     * @return minTotalDistance -> the shortest path length between two terms, -1 if no common ancestor is found
     * Computes the shortest path length (number of edges) between two terms
     * by computing their ancestor maps (using parent links only) and then finding a common
     * ancestor with the minimal total distance.
     * Returns -1 if no common ancestor is found.
     */
    private int computeShortestPath(GOTerm term1, GOTerm term2) {
        HashMap<GOTerm, Integer> ancestors1 = computeAncestors(term1);
        HashMap<GOTerm, Integer> ancestors2 = computeAncestors(term2);
        int minTotalDistance = Integer.MAX_VALUE;
        for (Map.Entry<GOTerm, Integer> entry : ancestors1.entrySet()) {
            GOTerm ancestor = entry.getKey();
            if (ancestors2.containsKey(ancestor)) {
                int totalDistance = entry.getValue() + ancestors2.get(ancestor);
                minTotalDistance = Math.min(minTotalDistance, totalDistance);
            }
        }
        return (minTotalDistance == Integer.MAX_VALUE) ? -1 : minTotalDistance;
    }

    private boolean isRelative(GOTerm term1, GOTerm term2) {
        // Relative if one is an ancestor of the other
        HashMap<GOTerm, Integer> ancestors1 = computeAncestors(term1);
        if (ancestors1.containsKey(term2)) {
            return true;
        }
        HashMap<GOTerm, Integer> ancestors2 = computeAncestors(term2);
        return ancestors2.containsKey(term1);
    }


    public static final class Builder {
        private RootType rootType;
        private int minSize;
        private int maxSize;
        private String output;
        private String overlapOut;

        public Builder() {
        }

        public Builder setRootType(RootType rootType) {
            this.rootType = rootType;
            return this;
        }

        public Builder setMinSize(int minSize) {
            this.minSize = minSize;
            return this;
        }

        public Builder setMaxSize(int maxSize) {
            this.maxSize = maxSize;
            return this;
        }

        public Builder setOverlapOut(String overlapOut) {
            this.overlapOut = overlapOut;
            return this;
        }

        public Builder setOutput(String output) {
            this.output = output;
            return this;
        }

        public GeneSetEnrichmentAnalysis build() {
            return new GeneSetEnrichmentAnalysis(this);
        }
    }

    public record DifferentialExpressionRecord(double fc, boolean signif) {
    }
}
