package gtf.treecollections;

import bamfeatures.SAMReadPair;
import gtf.GTFAnnotation;
import gtf.structs.Gene;
import gtf.types.StrandDirection;

import java.util.List;

public interface IntervalTreeForestManager {

    String currentChromosome = null;

    /**
     * @param chromosome the chromosome to be selected
     *                   Deletes the current tree and selects the tree corresponding to the chromosome, changes currentChromosome
     */
    void nextTree(String chromosome);

    /**
     * @param pair the pair to be checked
     * @return true if the current tree has a gene that is enclosed by the pair
     */
    boolean hasContainedGene(SAMReadPair pair);

    /**
     * @param pair the pair to be checked
     * @return the genes that are enclosing the pair (pair is inside the gene)
     */
    List<Gene> getGenesThatInclude(SAMReadPair pair);

    /**
     * @param pair the pair to be checked
     * @return The (minimum) distance to the nearest gene in the current tree
     */
    int getDistanceToNearestNeighborGene(SAMReadPair pair);

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic, MergedTranscriptomic or Intronic match in the current tree
     */
    boolean isAntisenseBetterThanAll(SAMReadPair pair);

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic or MergedTranscriptomic match in the current tree
     */
    boolean isAntisenseBetterThanIntronic(SAMReadPair pair);

    /**
     * @param pair the pair to be checked
     * @return true if the pair has a Transcriptomic match in the current tree
     */
    boolean isAntisenseBetterThanMergedTranscriptomic(SAMReadPair pair);

    /**
     * @param gtfAnnotation the annotation to be used
     *                      Initializes the IntervalTreeForestManager with the given annotation
     * @return
     */
    void init(GTFAnnotation gtfAnnotation);


}
