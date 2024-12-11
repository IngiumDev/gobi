package me.ingiumdev.gobi.gtf.structs;

import augmentedTree.IntervalTree;
import me.ingiumdev.gobi.gtf.types.StrandDirection;
import me.ingiumdev.gobi.readsimulator.Pair;

public class ChromosomeTreePair extends Pair<IntervalTree<Gene>, IntervalTree<Gene>> {
    public ChromosomeTreePair(IntervalTree<Gene> first, IntervalTree<Gene> second) {
        super(first, second);
    }

    public ChromosomeTreePair() {
        super(new IntervalTree<>(), new IntervalTree<>());
    }

    public IntervalTree<Gene> getPlusTree() {
        return this.getFirst();
    }

    public IntervalTree<Gene> getMinusTree() {
        return this.getSecond();
    }

    public IntervalTree<Gene> getTree(StrandDirection strand) {
        if (strand == StrandDirection.FORWARD) {
            return this.getPlusTree();
        } else if (strand == StrandDirection.REVERSE) {
            return this.getMinusTree();
        } else {
            throw new IllegalArgumentException("StrandDirection must be FORWARD or REVERSE");
        }
    }
}
