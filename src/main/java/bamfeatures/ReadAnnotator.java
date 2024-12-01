package bamfeatures;

import gtf.GTFAnnotation;
import gtf.structs.Gene;
import gtf.structs.Transcript;
import gtf.treecollections.*;
import gtf.types.StrandDirection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import parsers.GTFParser;
import readsimulator.Pair;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class ReadAnnotator {
    private final StrandDirection strandSpecificity;
    private SamReader samReader;
    private File outputFile;
    private IntervalTreeForestManager forestManager;
    private File gtfFile;
    private HashMap<String, SAMRecord> lookup;
    private List<SAMReadPair> readsToAnnotate;
    private String currentChromosome = "_";
    private PCRIndexManager pcrIndex;

    private ReadAnnotator(Builder builder) {
        samReader = builder.samReader;
        gtfFile = builder.gtfFile;
        outputFile = builder.outputFile;
        strandSpecificity = builder.strandSpecificity;
    }

    public void annotateReads() {
        // TODO: init forestManager
        Iterator<SAMRecord> it = samReader.iterator();
        if (strandSpecificity == StrandDirection.UNSPECIFIED) {
            forestManager = new StrandUnspecificForest();
            pcrIndex = new StrandUnSpecificPCRIndex();
        } else {
            //TODO: Possible migrate to tree pair instead of hasmap of strands
            forestManager = new StrandSpecificForest(strandSpecificity);
            pcrIndex = new StrandSpecificPCRIndex();
        }
        // TODO: better GTF file handling
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(String.valueOf(gtfFile));
        forestManager.init(gtfAnnotation);
        pcrIndex.initializePCRIndex();
        String referenceName;
        String readName;
        SAMRecord record;
        boolean isFirstOfPair;
        while (it.hasNext()) {
            record = it.next();
            referenceName = record.getReferenceName();
            if (!currentChromosome.equals(referenceName)) {
                processNewChromosome(referenceName);
            }
            if (isValidRead(record)) {
                readName = record.getReadName();
                if (lookup.containsKey(readName)) {
                    if (readName.equals("19773803")) {
                        System.out.println();

                    }
                    // Check if the read is the first or second of the pair
                    if (record.getFirstOfPairFlag()) {
                        readsToAnnotate.add(new SAMReadPair(record, lookup.get(readName)));
                        ReadAnnotation ra = new ReadAnnotation(readName);
                        ra.extractReadIntervals(record, lookup.get(readName));
                        ra.extractReadAlignmentStartEnd(record, lookup.get(readName));

                        ra.calculateClipping(record, lookup.get(readName));
                        ra.calculateMismatches(record, lookup.get(readName));
                        ra.calculateSplit();
                        //System.out.println(readName + "\tmm:" + ra.getMismatchCount() + "\tclipping:" + ra.getClippingSum() + "\tnsplit:" + ra.getSplitCount());
                        List<Gene> resultGenes = forestManager.getGenesThatInclude(ra);
                        if (!resultGenes.isEmpty()) {
                            // transcriptomic process
                            ra.setGenesThatInclude(resultGenes);
                            // transcriptomic process
                            if (!ra.areReadsConsistent()) {
//                                System.out.println(readName + "\tsplit-inconsistent:true");
                            } else {
                                ra.calculatePCRIndex(pcrIndex);

                                if (ra.findTranscriptomicMatches()) {
                                    // Gene-id,Gene-biotype:transcript-id1,transcript-id2|Gene-id,Gene-biotype:transcript-id1,transcript-id2

                                    StringBuilder sb = new StringBuilder();
                                    sb.append(readName + "\t" + "mm:" + ra.getMismatchCount() + "\tclipping:" + ra.getClippingSum() + "\tnsplit:" + ra.getSplitCount() + "\tgcount:" + ra.getGeneCount() + "\t");
                                    for (Pair<Gene, List<Transcript>> match : ra.getTranscriptomicMatches()) {
                                        sb.append(match.getFirst().getGeneID()).append(",").append(match.getFirst().getSource()).append(":");
                                        for (Transcript t : match.getSecond()) {
                                            sb.append(t.getTranscriptID()).append(",");
                                        }
                                        sb.deleteCharAt(sb.length() - 1);
                                        sb.append("|");
                                    }
                                    sb.deleteCharAt(sb.length() - 1);
                                    sb.append("\tpcrindex: ").append(ra.getPcrIndex());

                                    System.out.println(sb.toString());

                                } else {
                                    if (ra.findMergedTranscriptomicMatches()) {
//                                        System.out.print("");
                                    } else {

                                    }
                                }


                            }
                        } else {
                            // Check whether it contains a gene
                            if (!forestManager.hasContainedGene(ra)) {
                                // gdist
                                int gdist = forestManager.getDistanceToNearestNeighborGene(ra);
                                if (ra.areReadsConsistent()) {
                                    ra.calculatePCRIndex(pcrIndex);
//                                    System.out.println(readName + "\tmm:" + ra.getMismatchCount() + "\tclipping:" + ra.getClippingSum() + "\tnsplit:" + ra.getSplitCount() + "\tgcount:0" + "\tgdist:" + gdist+"\tpcrindex: "+ra.getPcrIndex());

                                } else {
//                                    System.out.println(readName + "\tsplit-inconsistent:true");
                                }
                            } else {
                                // skip
                            }
                        }


                    } else {
                        readsToAnnotate.add(new SAMReadPair(lookup.get(readName), record));
                        ReadAnnotation ra = new ReadAnnotation(readName);
                        ra.extractReadIntervals(lookup.get(readName), record);
                        ra.extractReadAlignmentStartEnd(lookup.get(readName), record);


                        ra.calculateClipping(lookup.get(readName), record);
                        ra.calculateMismatches(lookup.get(readName), record);
                        ra.calculateSplit();
//                            System.out.println(readName + "\tmm:" + ra.getMismatchCount() + "\tclipping:" + ra.getClippingSum() + "\tnsplit:" + ra.getSplitCount());
                        List<Gene> resultGenes = forestManager.getGenesThatInclude(ra);
                        if (!resultGenes.isEmpty()) {
                            ra.setGenesThatInclude(resultGenes);
                            // transcriptomic process
                            if (!ra.areReadsConsistent()) {
//                                System.out.println(readName + "\tsplit-inconsistent:true");
                            } else {
                                ra.calculatePCRIndex(pcrIndex);
                                if (ra.findTranscriptomicMatches()) {
                                    StringBuilder sb = new StringBuilder();
                                    sb.append(readName + "\t" + "mm:" + ra.getMismatchCount() + "\tclipping:" + ra.getClippingSum() + "\tnsplit:" + ra.getSplitCount() + "\tgcount:" + ra.getGeneCount() + "\t");
                                    for (Pair<Gene, List<Transcript>> match : ra.getTranscriptomicMatches()) {
                                        sb.append(match.getFirst().getGeneID()).append(",").append(match.getFirst().getSource()).append(":");
                                        for (Transcript t : match.getSecond()) {
                                            sb.append(t.getTranscriptID()).append(",");
                                        }
                                        sb.deleteCharAt(sb.length() - 1);
                                        sb.append("|");
                                    }
                                    sb.deleteCharAt(sb.length() - 1);
                                    sb.append("\tpcrindex: ").append(ra.getPcrIndex());

                                    System.out.println(sb.toString());
                                } else {
                                    if (ra.findMergedTranscriptomicMatches()) {
                                    } else {

                                    }
                                }

                            }
                        } else {
                            // Check whether it contains a gene
                            if (!forestManager.hasContainedGene(ra)) {
                                // gdist
                                int gdist = forestManager.getDistanceToNearestNeighborGene(ra);
                                if (ra.areReadsConsistent()) {
                                    ra.calculatePCRIndex(pcrIndex);
//                                    System.out.println(readName + "\tmm:" + ra.getMismatchCount() + "\tclipping:" + ra.getClippingSum() + "\tnsplit:" + ra.getSplitCount() + "\tgcount:0" + "\tgdist:" + gdist+"\tpcrindex: "+ra.getPcrIndex());

                                } else {
//                                    System.out.println(readName + "\tsplit-inconsistent:true");
                                }
                            } else {
                                // skip
                            }
                        }


                    }


                } else {
                    lookup.put(readName, record);
                }
            }


        }
    }

    public boolean areReadsSameStrand(SAMRecord record) {
        return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
    }

    public boolean isValidRead(SAMRecord record) {
        return !record.getReadUnmappedFlag() && record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !areReadsSameStrand(record) && areReadsSameChromosome(record) && !record.getSupplementaryAlignmentFlag();
    }

    public boolean areReadsSameChromosome(SAMRecord record) {
        return record.getReferenceName().equals(record.getMateReferenceName());
    }

    private void processNewChromosome(String referenceName) {
        lookup = new HashMap<>();
        currentChromosome = referenceName;
        forestManager.nextTree(referenceName);
        // TODO: implement readprocessing
        readsToAnnotate = new ArrayList<>();
        pcrIndex.nextChromosome();

    }


    public static final class Builder {
        private SamReader samReader;
        private File gtfFile;
        private File outputFile;
        private StrandDirection strandSpecificity;

        public Builder() {
        }

        public Builder setSamReader(SamReader val) {
            samReader = val;
            return this;
        }

        public Builder setGtfFile(File val) {
            gtfFile = val;
            return this;
        }

        public Builder setOutputFile(File val) {
            outputFile = val;
            return this;
        }

        public Builder setStrandSpecificity(StrandDirection val) {
            strandSpecificity = val;
            return this;
        }

        public ReadAnnotator build() {
            return new ReadAnnotator(this);
        }
    }
}
