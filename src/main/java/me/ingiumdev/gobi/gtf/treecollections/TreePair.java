package me.ingiumdev.gobi.gtf.treecollections;

import augmentedTree.IntervalTree;
import me.ingiumdev.gobi.gtf.structs.Gene;
import me.ingiumdev.gobi.readsimulator.IdenticalPair;

public class TreePair extends IdenticalPair<IntervalTree<Gene>> {


    public TreePair(IntervalTree<Gene> first, IntervalTree<Gene> second) {
        super(first, second);
    }

    public TreePair() {
        super(new IntervalTree<>(), new IntervalTree<>());
    }
}
