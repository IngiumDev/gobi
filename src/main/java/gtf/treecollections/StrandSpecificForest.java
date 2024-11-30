package gtf.treecollections;

import augmentedTree.IntervalTree;
import bamfeatures.SAMReadPair;
import gtf.GTFAnnotation;
import gtf.structs.Gene;
import gtf.types.StrandDirection;

import java.util.*;

public class StrandSpecificForest implements IntervalTreeForestManager{
    Map<String, TreePair> chromosomeToGeneTree;
    StrandDirection strandSpecificity;
    TreePair currentTreePair;
    List<Gene> resultGenes;

    public StrandSpecificForest(StrandDirection strandSpecificity) {
        this.strandSpecificity = strandSpecificity;
    }

    /**
     * @param chromosome the chromosome to be selected
     *                   Deletes the current tree and selects the tree corresponding to the chromosome, changes currentChromosome
     */
    @Override
    public void nextTree(String chromosome) {
        currentTreePair = chromosomeToGeneTree.get(chromosome);
    }

    /**
     * @param pair the pair to be checked
     * @return true if the current tree has a gene that is enclosed by the pair
     */
    @Override
    public boolean hasContainedGene(SAMReadPair pair) {
        resultGenes = new ArrayList<>();
        // TODO: NEED TO CHECK THE RIGHT ONE, stranded needs to get the correct tree depending on getnegative of the first pair
        // TODO refactor to use different pair doesn't work for turned reads
        int minstart = Math.min(pair.getFirst().getAlignmentStart(), pair.getSecond().getAlignmentStart());
        int maxend = Math.max(pair.getFirst().getAlignmentEnd(), pair.getSecond().getAlignmentEnd());
        currentTreePair.getFirst().getIntervalsSpannedBy(minstart, maxend, resultGenes);
        return !resultGenes.isEmpty();
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
        chromosomeToGeneTree = new HashMap<>();
        // TODO: Biotype extracting
        for (Gene gene : gtfAnnotation.getGenes().values()) {
            String chromosome = gene.getSeqname();
            StrandDirection strand = gene.getStrand();
            if (!chromosomeToGeneTree.containsKey(chromosome)) {
                chromosomeToGeneTree.put(chromosome, new TreePair());
            }
            if (strandSpecificity == StrandDirection.FORWARD) {
                if (strand == StrandDirection.FORWARD) {
                    chromosomeToGeneTree.get(chromosome).getFirst().add(gene);
                } else {
                    // StrandDirection.REVERSE
                    chromosomeToGeneTree.get(chromosome).getSecond().add(gene);
                }
            } else {
                // StrandDirection.REVERSE
                if (strand == StrandDirection.REVERSE) {
                    chromosomeToGeneTree.get(chromosome).getFirst().add(gene);
                } else {
                    // StrandDirection.FORWARD
                    chromosomeToGeneTree.get(chromosome).getSecond().add(gene);
                }
            }
        }

    }
}
