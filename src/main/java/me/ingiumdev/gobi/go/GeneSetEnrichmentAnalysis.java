package me.ingiumdev.gobi.go;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class GeneSetEnrichmentAnalysis {
    int numDiffExpressedRootGenes;
    int numDiffExpressedSigRootGenes;
    private DAG graph;
    private Mapping mapping;
    private RootType rootType;
    private int minSize;
    private int maxSize;
    private Map<String, DifferentialExpressionRecord> differentialExpressionInput;
    private Set<Integer> SOTterms;

    public GeneSetEnrichmentAnalysis(RootType rootType) {
    }

    private GeneSetEnrichmentAnalysis(Builder builder) {
        rootType = builder.rootType;
        minSize = builder.minSize;
        maxSize = builder.maxSize;
    }

    public void initMapping(String mappingPath) {
        mapping = Mapping.createEnsemblMapping(mappingPath);
    }

    public void initDAG(String oboPath) {
        graph = new DAG(oboPath, mapping, rootType);
    }

    public void initDifferentialExpression(String path) {
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
    }

    public List<GOAnalysisEntry> performEnrichment() {
        calcNumDiffExpressedRootGenes(graph, differentialExpressionInput);
        differentialExpressionInput.size();
        // Replace with the actual transformation logic
        List<GOAnalysisEntry> analysisEntries = graph.getEntries().values().parallelStream()
                .filter(goTerm -> minSize <= goTerm.getAssociatedGenes().size()
                        && goTerm.getAssociatedGenes().size() <= maxSize)
                .map(this::analyzeGOEntry).sorted(Comparator.comparingInt(GOAnalysisEntry::getId)).collect(Collectors.toList());
        // fdr
        // output the results and order by id ascending
        //term	name	size	is_true	noverlap	hg_pval	hg_fdr	fej_pval	fej_fdr	ks_stat	ks_pval	ks_fdr	shortest_path_to_a_true
        System.out.println("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true");
        analysisEntries.forEach(System.out::println);

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

// N is total genes in root interescted with the genes in the differential expression
        //K is the number of genes in the differential expression
        // n = size of the GO term
        // k =numoverlap
        return entry;
    }

    public static final class Builder {
        private RootType rootType;
        private int minSize;
        private int maxSize;

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

        public GeneSetEnrichmentAnalysis build() {
            return new GeneSetEnrichmentAnalysis(this);
        }
    }

    record DifferentialExpressionRecord(double fc, boolean signif) {
    }
}
