package me.ingiumdev.gobi.go;

import me.ingiumdev.gobi.parsers.OboParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

public class DAG {
    private static final Logger log = LoggerFactory.getLogger(DAG.class);
    private final HashMap<Integer, GOTerm> entries;
    private GOTerm root;
    private List<GOTerm> leaves;

    public DAG(String path, Mapping mapping, RootType rootType) {
        OboParser parser = new OboParser(path, rootType);
        long start = System.currentTimeMillis();
        entries = parser.parse();
        // loaded entries
        log.info("Loaded {} entries in {}ms", entries.size(), System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        resolveParents();
        log.info("Resolved parents in {}ms", System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        resolveLeaves(); // root should have 15497 genes in ensembl
        log.info("Resolved leaves in {}ms", System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        resolveGenes(mapping);
        log.info("Resolved genes in {}ms", System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        propagateGenes(root);
        log.info("Propagated genes in {}ms", System.currentTimeMillis() - start);
        log.info("Root contains {} genes", root.getAssociatedGenes().size());
    }

    public GOTerm getRoot() {
        return root;
    }

    public HashMap<Integer, GOTerm> getEntries() {
        return entries;
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

    public List<GOTerm> getLeaves() {
        return leaves;
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
