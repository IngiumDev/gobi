package bamfeatures;

import gtf.GTFAnnotation;
import gtf.structs.Gene;
import gtf.treecollections.*;
import gtf.types.StrandDirection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import parsers.GTFParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;

public class ReadAnnotator {
    private final StrandDirection strandSpecificity;
    private final SamReader samReader;
    private final File outputFile;
    private final File gtfFile;
    private final List<ReadAnnotation> allReadAnnotations = new ArrayList<>();
    private IntervalTreeForestManager forestManager;
    private HashMap<String, SAMRecord> lookup;
    private List<SAMReadPair> readsToAnnotate;
    private String currentChromosome = "_";
    private PCRIndexManager pcrIndex;
    private BufferedWriter writer;
    private boolean returnAll = false;

    private ReadAnnotator(Builder builder) {
        samReader = builder.samReader;
        gtfFile = builder.gtfFile;
        outputFile = builder.outputFile;
        strandSpecificity = builder.strandSpecificity;
    }

    public List<ReadAnnotation> annotateAndReturnReads() {

        returnAll = true;
        annotateReads();
        return allReadAnnotations;
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
        if (outputFile != null) {
            try {
                writer = new BufferedWriter(new FileWriter(outputFile));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

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
                    // Check if the read is the first or second of the pair
                    if (record.getFirstOfPairFlag()) {
                        readsToAnnotate.add(new SAMReadPair(record, lookup.get(readName)));
                    } else {
                        readsToAnnotate.add(new SAMReadPair(lookup.get(readName), record));
                    }
                } else {
                    lookup.put(readName, record);
                }
            }
        }
        processLastChromosome();

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
        // TODO: implement readprocessing
        if (readsToAnnotate != null) {
            processAndWriteReads();

        }
        readsToAnnotate = new ArrayList<>();
        lookup = new HashMap<>();
        currentChromosome = referenceName;
        forestManager.nextTree(referenceName);

        pcrIndex.nextChromosome();

    }

    private void outputReads(List<ReadAnnotation> readAnnotations, File outputFile) {

        for (ReadAnnotation readAnnotation : readAnnotations) {
            if (readAnnotation != null) {
                try {
                    writer.write(readAnnotation.output());
                    writer.newLine();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }

        }

    }

    private void processAndWriteReads() {
        if (outputFile != null) {
            AtomicBoolean producerDone = new AtomicBoolean(false);
            ConcurrentLinkedQueue<ReadAnnotation> queue = new ConcurrentLinkedQueue<>();
            Thread producerThread = new Thread(() -> {
                for (SAMReadPair samReadPair : readsToAnnotate) {
                    ReadAnnotation readAnnotation = processRead(samReadPair);
                    if (readAnnotation != null) {
                        queue.add(readAnnotation);
                        if (returnAll) {
                            allReadAnnotations.add(readAnnotation);
                        }
                    }
                }
                producerDone.set(true);
            });
            Thread consumerThread = new Thread(() -> {
                try {
                    while (!producerDone.get() || !queue.isEmpty()) {
                        ReadAnnotation data = queue.poll();
                        if (data != null) {
                            writer.write(data.output());
                            writer.newLine();
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
            producerThread.start();
            consumerThread.start();
            try {
                producerThread.join();
                consumerThread.join();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }

        } else if (returnAll) {
            for (SAMReadPair samReadPair : readsToAnnotate) {
                ReadAnnotation readAnnotation = processRead(samReadPair);
                if (readAnnotation != null) {
                    allReadAnnotations.add(readAnnotation);
                }
            }
        }
    }

    private ReadAnnotation processRead(SAMReadPair samReadPair) {
        SAMRecord first = samReadPair.getFirst();
        SAMRecord second = samReadPair.getSecond();
        ReadAnnotation readAnnotation = new ReadAnnotation(first.getReadName());
        readAnnotation.extractReadAlignmentStartEnd(first, second);
        List<Gene> resultGenes = forestManager.getGenesThatInclude(readAnnotation);
        if (!resultGenes.isEmpty()) {
            readAnnotation.setGenesThatInclude(resultGenes);
            readAnnotation.extractReadIntervals(first, second);
            if (readAnnotation.areReadsConsistent()) {
                calculateBasicReadInfo(readAnnotation, first, second);
                // Specifics
                if (!readAnnotation.findTranscriptomicMatches()) {
                    readAnnotation.findMergedTranscriptomicMatches();
                }
                return readAnnotation;
            } else {
                return readAnnotation;
            }
        } else if (!forestManager.hasContainedGene(readAnnotation)) {
            // Check whether it contains a gene
            readAnnotation.extractReadIntervals(first, second);
            if (readAnnotation.areReadsConsistent()) {
                calculateBasicReadInfo(readAnnotation, first, second);
                // Specifics
                readAnnotation.calculateGeneDistance(forestManager);
                readAnnotation.processAntisense(forestManager);
                return readAnnotation;
            } else {
                return readAnnotation;
            }

        }
        return null; // if contained a gene and was not included
    }

    private void calculateBasicReadInfo(ReadAnnotation readAnnotation, SAMRecord first, SAMRecord second) {
        readAnnotation.calculateClipping(first, second);
        readAnnotation.calculateMismatches(first, second);
        readAnnotation.calculateSplit();
        readAnnotation.calculatePCRIndex(pcrIndex);
    }


    private void processLastChromosome() {
        processAndWriteReads();
        try {
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
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
