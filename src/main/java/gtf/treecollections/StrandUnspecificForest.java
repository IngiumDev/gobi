package gtf.treecollections;

import augmentedTree.IntervalTree;
import bamfeatures.SAMReadPair;
import gtf.GTFAnnotation;
import gtf.structs.Gene;

import java.util.HashMap;
import java.util.List;

public class StrandUnspecificForest implements IntervalTreeForestManager {
    HashMap<String, IntervalTree<Gene>> chromosomeToGeneTree = new HashMap<>();

    List<Gene> resultGenes;

    /**
     * @param chromosome the chromosome to be selected
     *                   Deletes the current tree and selects the tree corresponding to the chromosome, changes currentChromosome
     */
    @Override
    public void nextTree(String chromosome) {

    }


    /**
     * @param pair the pair to be checked
     * @return true if the current tree has a gene that is enclosed by the pair
     */
    @Override
    public boolean hasContainedGene(SAMReadPair pair) {
        return false;
    }

    /**
     * @param pair the pair to be checked
     * @return the genes that are enclosing the pair (pair is inside the gene)
     */
    @Override
    public List<Gene> getGenesThatInclude(SAMReadPair pair) {
        return List.of();
    }

    /**
     * @param pair the pair to be checked
     * @return The (minimum) distance to the nearest gene in the current tree
     */
    @Override
    public int getDistanceToNearestNeighborGene(SAMReadPair pair) {
        return 0;
    }

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic, MergedTranscriptomic or Intronic match in the current tree
     */
    @Override
    public boolean isAntisenseBetterThanAll(SAMReadPair pair) {
        return false;
    }

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic or MergedTranscriptomic match in the current tree
     */
    @Override
    public boolean isAntisenseBetterThanIntronic(SAMReadPair pair) {
        return false;
    }

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic match in the current tree
     */
    @Override
    public boolean isAntisenseBetterThanMergedTranscriptomic(SAMReadPair pair) {
        return false;
    }

    /**
     * @param gtfAnnotation the annotation to be used
     *                      Initializes the IntervalTreeForestManager with the given annotation
     * @return
     */
    @Override
    public void init(GTFAnnotation gtfAnnotation) {
    }

}
