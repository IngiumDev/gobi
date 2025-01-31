package me.ingiumdev.gobi.go;

import me.ingiumdev.gobi.parsers.OboParser;

import java.util.*;

public class DAG {
    private final HashMap<Integer, GOTerm> entries;
    private GOTerm root;
    private List<GOTerm> leaves;

    public DAG(String path, Mapping mapping, RootType rootType) {
        OboParser parser = new OboParser(path, rootType);
        entries = parser.parse();
        resolveParents();
        resolveLeaves(); // root should have 15497 genes in ensembl
        resolveGenes(mapping);
        propagateGenes(root);
        System.out.println();
    }

    public GOTerm getRoot() {
        return root;
    }

    public HashMap<Integer, GOTerm> getEntries() {
        return entries;
    }

    public static void main(String[] args) {
        long start = System.currentTimeMillis();
        Mapping mapping = Mapping.createEnsemblMapping("~/IdeaProjects/gobi/data/GOEnrich/goa_human_ensembl.tsv");
        DAG dag = new DAG("~/IdeaProjects/gobi/data/GOEnrich/go.obo", mapping, RootType.BIOLOGICAL_PROCESS);
        System.out.println(System.currentTimeMillis() - start);
    }

    private void propagateGenes(GOTerm node) {
        for (GOTerm child : node.getChildren()) {
            propagateGenes(child);
            node.getAssociatedGenes().addAll(child.getAssociatedGenes());
        }
    }

    public GOTerm getTerm(int id) {
        return entries.get(id);
    }
    private void propagateUpward(GOTerm leaf) {
        List<GOTerm> parents = leaf.getParents();
        for (GOTerm parent : parents) {
            parent.getAssociatedGenes().addAll(leaf.getAssociatedGenes());
            propagateUpward(parent);
        }
    }

    private void resolveGenes(Mapping mapping) {
        for (Map.Entry<String, Set<Integer>> entry : mapping.getMap().entrySet()) {
            String gene = entry.getKey();
            Set<Integer> associatedTerms = entry.getValue();

            for (int term : associatedTerms) {
                GOTerm goTerm = entries.get(term);
                if (goTerm != null) {
                    goTerm.addGene(gene);
                }
            }
        }
    }

    public void resolveParents() {
        GOTerm parent;
        for (GOTerm entry : entries.values()) {
            for (int parentID : entry.getTempParentIDs()) {
                parent = entries.get(parentID);
                if (parent != null) {
                    entry.addParent(parent);
                    parent.addChild(entry);
                }
            }
            if (entry.getTempParentIDs().isEmpty()) {
                root = entry;
            }
            entry.clearTempParentIDs();
        }
    }

    public void resolveLeaves() {
        leaves = new ArrayList<>();
        for (GOTerm entry : entries.values()) {
            if (entry.getChildren().isEmpty()) {
                leaves.add(entry);
            }
        }
    }


}
