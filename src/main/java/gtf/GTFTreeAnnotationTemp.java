package gtf;

import augmentedTree.IntervalTree;
import gtf.structs.ChromosomeTreePair;
import gtf.structs.Gene;
import gtf.types.StrandDirection;

import java.util.HashMap;

public class GTFTreeAnnotationTemp {
    // TODO Refactor GTF Parsing, make this calculate when needed only
    HashMap<String, ChromosomeTreePair> chromosomeTreePairHashMap = new HashMap<>();

    public GTFTreeAnnotationTemp(GTFAnnotation gtfAnnotation) {
        long time = System.currentTimeMillis();
        for (Gene gene : gtfAnnotation.getGenes().values()) {
            String chromosome = gene.getSeqname();
            if (!chromosomeTreePairHashMap.containsKey(chromosome)) {
                chromosomeTreePairHashMap.put(chromosome, new ChromosomeTreePair());

            }
            ChromosomeTreePair chromosomeTreePair = chromosomeTreePairHashMap.get(chromosome);
            if (gene.getStrand() == StrandDirection.FORWARD) {
                chromosomeTreePair.getPlusTree().add(gene);
            } else if (gene.getStrand() == StrandDirection.REVERSE) {
                chromosomeTreePair.getMinusTree().add(gene);
            }
        }
        System.out.println("LOG: Total time to build chromosome tree pair: " + (System.currentTimeMillis() - time) + " ms");
    }

    public IntervalTree<Gene> getTree(String chromosome, StrandDirection strand) {
        if (chromosomeTreePairHashMap.containsKey(chromosome)) {
            return chromosomeTreePairHashMap.get(chromosome).getTree(strand);
        }
        return null;
    }
}
