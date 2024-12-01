package gtf.treecollections;

import bamfeatures.ReadAnnotation;
import gtf.GTFAnnotation;
import gtf.structs.Gene;
import gtf.types.StrandDirection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class StrandSpecificForest implements IntervalTreeForestManager {
    Map<String, TreePair> chromosomeToGeneTree;
    StrandDirection strandSpecificity;
    TreePair currentTreePair;

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
        // TODO: delete the current tree
    }

    /**
     * @param pair the pair to be checked
     * @return true if the current tree has a gene that is enclosed by the pair
     */
    @Override
    public boolean hasContainedGene(ReadAnnotation pair) {
        List<Gene> resultGenes = new ArrayList<>();
        if (!pair.isReadStrandNegative()) {
            currentTreePair.getFirst().getIntervalsSpannedBy(pair.getAlignmentStart(), pair.getAlignmentEnd(), resultGenes);
        } else {
            currentTreePair.getSecond().getIntervalsSpannedBy(pair.getAlignmentStart(), pair.getAlignmentEnd(), resultGenes);
        }
        // No need to return the list, just check if it is empty
        return !resultGenes.isEmpty();
    }

    /**
     * @param pair the pair to be checked
     * @return the genes that are enclosing the pair (pair is inside the gene)
     */
    @Override
    public List<Gene> getGenesThatInclude(ReadAnnotation pair) {
        List<Gene> resultGenes = new ArrayList<>();
        if (!pair.isReadStrandNegative()) {
            currentTreePair.getFirst().getIntervalsSpanning(pair.getAlignmentStart(), pair.getAlignmentEnd(),resultGenes);
        } else {
            currentTreePair.getSecond().getIntervalsSpanning(pair.getAlignmentStart(), pair.getAlignmentEnd(),resultGenes);
        }
        // WARNING: Maybe size() == 0;
        return resultGenes;
    }

    /**
     * @param pair the pair to be checked
     * @return The (minimum) distance to the nearest gene in the current tree
     */
    @Override
    public int getDistanceToNearestNeighborGene(ReadAnnotation pair) {
        int leftDistance;
        int rightDistance;
        if (!pair.isReadStrandNegative()) {
            leftDistance = IntervalTreeForestManager.getLeftDistance(currentTreePair.getFirst(), pair);
            rightDistance = IntervalTreeForestManager.getRightDistance(currentTreePair.getFirst(), pair);
        } else {
            leftDistance = IntervalTreeForestManager.getLeftDistance(currentTreePair.getSecond(), pair);
            rightDistance = IntervalTreeForestManager.getRightDistance(currentTreePair.getSecond(), pair);
        }
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
