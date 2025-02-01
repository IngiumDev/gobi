package me.ingiumdev.gobi.go;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class GeneSetEnrichmentAnalysis {
    private static final Logger log = LoggerFactory.getLogger(GeneSetEnrichmentAnalysis.class);
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

    /**
     * Applies the Benjamini–Hochberg correction on a list of entries.
     *
     * @param entries      The list of entries to adjust.
     * @param pValueGetter A function that extracts the raw p-value from an entry.
     * @param fdrSetter    A consumer that sets the corrected FDR value on an entry.
     * @param <T>          The type of the entries (e.g., GOAnalysisEntry).
     */
    public static <T> void applyBenjaminiHochbergCorrection(List<T> entries, Function<T, Double> pValueGetter, BiConsumer<T, Double> fdrSetter) {
        int m = entries.size();

        // Sort the entries in ascending order based on the extracted p-value.
        List<T> sortedEntries = entries.stream().sorted(Comparator.comparingDouble(pValueGetter::apply)).collect(Collectors.toList());

        double[] adjustedFDR = new double[m];

        // Iterate backwards to compute the BH-adjusted p-values.
        for (int i = m - 1; i >= 0; i--) {
            double pValue = pValueGetter.apply(sortedEntries.get(i));
            int rank = i + 1; // Ranks are 1-based
            double bhValue = pValue * m / rank;

            if (i == m - 1) {
                adjustedFDR[i] = bhValue;
            } else {
                // Ensure monotonicity by taking the minimum with the previously computed value.
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
        // Replace with the actual transformation logic
        List<GOAnalysisEntry> analysisEntries = graph.getEntries().values().parallelStream().filter(goTerm -> minSize <= goTerm.getAssociatedGenes().size() && goTerm.getAssociatedGenes().size() <= maxSize).map(this::analyzeGOEntry).sorted(Comparator.comparingInt(GOAnalysisEntry::getId)).collect(Collectors.toList());
        // calciulated basic encirhcment statistics in
        log.info("Calculated basic enrichment statistics in {} ms", System.currentTimeMillis() - start);
        // Do BH for each p-value
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

    // This method should be implemented to create a GOAnalysisEntry from a GOTerm
    private GOAnalysisEntry analyzeGOEntry(GOTerm goTerm) {
        // Process the GO term to generate a GOAnalysisEntry
        GOAnalysisEntry entry = new GOAnalysisEntry(goTerm);
        // calculate size
        entry.calculateSize(graph, differentialExpressionInput);
        entry.calculatehg_pvalue(numDiffExpressedRootGenes, numDiffExpressedSigRootGenes, goTerm);
        entry.calculatefej_pvalue(numDiffExpressedRootGenes, numDiffExpressedSigRootGenes, goTerm);
        entry.calculateKS(graph, differentialExpressionInput);
        entry.calculateShortestPathToTrue(SOTterms, graph, goTerm);
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

    /*
    • [-overlapout overlap_out_tsv]: optional parameter that specifies an output file.
    if set, information about DAG entries, with shared mapped genes is written into this
    file. The [overlapout tsv file must have the following columns:
    – term1: GO-id (example: GO:1902554) of the first of the two overlapping DAG
    entries
    – term2: GO-id of the second of the two overlapping DAG entries
    – is_relative: true if the associated DAG entry to term1 is ascendant or descendent of the one associated to term2, false otherwise.
    – path_length: the length of shortest path between the two DAG entries. The
    length is defined as the minimal number of edges between term1 and term2.
    Hint: there may exist a shorter path between relatives than the direct one.
    – num_overlapping: the number of gene ids associated to both DAG entries
    4
    – max_ov_percent: the maximum percentage (a float value between 0.0 and 100.0)
    of the shared gene ids to all associated gene ids to term1 or term2
    Hint: output only the GO entry pairs both passing the minsize, maxsize criteria.
    You find an example output for overlapout in go_bp_mapping_go_50_500.overlapout
    for the parameters:
     */
    public void performOverlapAnalysis() {
        long start = System.currentTimeMillis();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(overlapOut))) {
            writer.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent");
            // Retrieve and filter GOTerm objects based on the number of associated genes.
            List<GOTerm> filteredTerms = graph.getEntries().values().stream()
                    .filter(goTerm -> minSize <= goTerm.getAssociatedGenes().size() &&
                            goTerm.getAssociatedGenes().size() <= maxSize)
                    .toList();

            // Process unique pairs in parallel using an index-based parallel stream.
            IntStream.range(0, filteredTerms.size()).parallel().forEach(i -> {
                for (int j = i + 1; j < filteredTerms.size(); j++) {
                    // Check if there is an overlap between the two GO terms.
                    int numSharedGenes = 0;
                    for (String gene : filteredTerms.get(i).getAssociatedGenes()) {
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

    private String processOverlapPair(GOTerm goTerm1, GOTerm goTerm2, int numSharedGenes) {
        boolean isRelative = false;
        int pathLength = 0;

        double maxOvPercent = (double) numSharedGenes / Math.min(goTerm1.getAssociatedGenes().size(), goTerm2.getAssociatedGenes().size()) * 100;
        return "GO:" + goTerm1.getFullID() + "\t" + "GO:" + goTerm2.getFullID() + "\t" + isRelative + "\t" + pathLength + "\t" + numSharedGenes + "\t" + maxOvPercent;
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

    record DifferentialExpressionRecord(double fc, boolean signif) {
    }
}
