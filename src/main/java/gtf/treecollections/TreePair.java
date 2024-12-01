package gtf.treecollections;

import augmentedTree.IntervalTree;
import gtf.structs.Gene;
import readsimulator.IdenticalPair;

public class TreePair extends IdenticalPair<IntervalTree<Gene>> {


    public TreePair(IntervalTree<Gene> first, IntervalTree<Gene> second) {
        super(first, second);
    }

    public TreePair() {
        super(new IntervalTree<>(), new IntervalTree<>());
    }
}
