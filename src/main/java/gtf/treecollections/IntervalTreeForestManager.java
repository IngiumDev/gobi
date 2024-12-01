package gtf.treecollections;

import augmentedTree.IntervalTree;
import bamfeatures.ReadAnnotation;
import gtf.GTFAnnotation;
import gtf.structs.Gene;

import java.util.ArrayList;
import java.util.List;

public interface IntervalTreeForestManager {

    String currentChromosome = null;

    static int getLeftDistance(IntervalTree<Gene> tree, ReadAnnotation pair) {
        List<Gene> genes = new ArrayList<>();
        tree.getIntervalsLeftNeighbor(pair.getAlignmentStart(), pair.getAlignmentEnd(), genes);
        if (!genes.isEmpty()) {
            Gene nearestGene = genes.getFirst();
            if (nearestGene.getStop() >= pair.getAlignmentStart()) {
                return 0;
            } else {
                return pair.getAlignmentStart() - nearestGene.getStop() - 1;
            }
        }
        return Integer.MAX_VALUE;
    }

    static int getRightDistance(IntervalTree<Gene> tree, ReadAnnotation pair) {
        List<Gene> genes = new ArrayList<>();
        tree.getIntervalsRightNeighbor(pair.getAlignmentStart(), pair.getAlignmentEnd(), genes);
        if (!genes.isEmpty()) {
            Gene nearestGene = genes.getFirst();
            if (nearestGene.getStart() <= pair.getAlignmentEnd()) {
                return 0;
            } else {
                return nearestGene.getStart() - pair.getAlignmentEnd() - 1;
            }
        }
        return Integer.MAX_VALUE;
    }

    /**
     * @param chromosome the chromosome to be selected
     *                   Deletes the current tree and selects the tree corresponding to the chromosome, changes currentChromosome
     */
    void nextTree(String chromosome);

    /**
     * @param pair the pair to be checked
     * @return true if the current tree has a gene that is enclosed by the pair
     */
    boolean hasContainedGene(ReadAnnotation pair);

    /**
     * @param pair the pair to be checked
     * @return the genes that are enclosing the pair (pair is inside the gene)
     */
    List<Gene> getGenesThatInclude(ReadAnnotation pair);

    /**
     * @param pair the pair to be checked
     * @return The (minimum) distance to the nearest gene in the current tree
     */
    int getDistanceToNearestNeighborGene(ReadAnnotation pair);

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic, MergedTranscriptomic or Intronic match in the current tree
     */
    boolean isAntisenseBetterThanAll(ReadAnnotation pair);

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic or MergedTranscriptomic match in the current tree
     */
    boolean isAntisenseBetterThanIntronic(ReadAnnotation pair);

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic match in the current tree
     */
    boolean isAntisenseBetterThanMergedTranscriptomic(ReadAnnotation pair);

    /**
     * @param gtfAnnotation the annotation to be used
     *                      Initializes the IntervalTreeForestManager with the given annotation
     * @return
     */
    void init(GTFAnnotation gtfAnnotation);


}
