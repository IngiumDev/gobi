package gtf.treecollections;

import augmentedTree.IntervalTree;
import bamfeatures.ReadAnnotation;
import gtf.GTFAnnotation;
import gtf.structs.Gene;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class StrandUnspecificForest implements IntervalTreeForestManager {
    HashMap<String, IntervalTree<Gene>> chromosomeToGeneTree = new HashMap<>();
    IntervalTree<Gene> currentTree;
    List<Gene> resultGenes;

    /**
     * @param chromosome the chromosome to be selected
     *                   Deletes the current tree and selects the tree corresponding to the chromosome, changes currentChromosome
     */
    @Override
    public void nextTree(String chromosome) {
        currentTree = chromosomeToGeneTree.get(chromosome);
        // TODO: delete the current tree
    }


    /**
     * @param pair the pair to be checked
     * @return true if the current tree has a gene that is enclosed by the pair
     */
    @Override
    public boolean hasContainedGene(ReadAnnotation pair) {
        List<Gene> resultGenes = new ArrayList<>();
        currentTree.getIntervalsSpannedBy(pair.getAlignmentStart(), pair.getAlignmentEnd(), resultGenes);

        return !resultGenes.isEmpty();
    }

    /**
     * @param pair the pair to be checked
     * @return the genes that are enclosing the pair (pair is inside the gene)
     */
    @Override
    public List<Gene> getGenesThatInclude(ReadAnnotation pair) {
        List<Gene> resultGenes = new ArrayList<>();
        currentTree.getIntervalsSpanning(pair.getAlignmentStart(), pair.getAlignmentEnd(), resultGenes);
        return resultGenes;
    }

    /**
     * @param pair the pair to be checked
     * @return The (minimum) distance to the nearest gene in the current tree
     */
    @Override
    public int getDistanceToNearestNeighborGene(ReadAnnotation pair) {
        int leftDistance = IntervalTreeForestManager.getLeftDistance(currentTree, pair);
        int rightDistance = IntervalTreeForestManager.getRightDistance(currentTree, pair);
        return Math.min(leftDistance, rightDistance);
    }

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic, MergedTranscriptomic or Intronic match in the current tree
     */
    @Override
    public boolean isAntisenseBetterThanAll(ReadAnnotation pair) {
        return false;
    }

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic or MergedTranscriptomic match in the current tree
     */
    @Override
    public boolean isAntisenseBetterThanIntronic(ReadAnnotation pair) {
        return false;
    }

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic match in the current tree
     */
    @Override
    public boolean isAntisenseBetterThanMergedTranscriptomic(ReadAnnotation pair) {
        return false;
    }

    /**
     * @param gtfAnnotation the annotation to be used
     *                      Initializes the IntervalTreeForestManager with the given annotation
     * @return
     */
    @Override
    public void init(GTFAnnotation gtfAnnotation) {
        chromosomeToGeneTree = new HashMap<>();
        for (Gene gene : gtfAnnotation.getGenes().values()) {
            if (!chromosomeToGeneTree.containsKey(gene.getSeqname())) {
                chromosomeToGeneTree.put(gene.getSeqname(), new IntervalTree<>());
            }
            chromosomeToGeneTree.get(gene.getSeqname()).add(gene);
        }
    }

}
