package bamfeatures;

import gtf.GTFAnnotation;
import gtf.structs.Interval;
import gtf.treecollections.IntervalTreeForestManager;
import gtf.treecollections.StrandSpecificForest;
import gtf.treecollections.StrandUnspecificForest;
import gtf.types.StrandDirection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import parsers.GTFParser;

import java.io.File;
import java.util.*;

public class ReadAnnotator {
    private SamReader samReader;
    private File outputFile;
    private final StrandDirection strandSpecificity;
    private IntervalTreeForestManager forestManager;
    private File gtfFile;
    private HashMap<String, SAMRecord> lookup;
    private List<SAMReadPair> readsToAnnotate;
    private String currentChromosome = "_";
    private HashMap<Set<Interval>, Integer> pcrIndex;

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
        } else {
            //TODO: Possible migrate to tree pair instead of hasmap of strands
            forestManager = new StrandSpecificForest(strandSpecificity);
        }
        // TODO: better GTF file handling
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(String.valueOf(gtfFile));
        forestManager.init(gtfAnnotation);
        System.out.println();
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
                        ReadAnnotation ra = new ReadAnnotation();
                        ra.extractReadIntervals(record, lookup.get(readName));
                        ra.extractReadAlignmentStartEnd(record, lookup.get(readName));
                        if (readName.equals("19216203")) {
                            System.out.println();
                        }
                        if (!ra.areReadsConsistent()) {
                            System.out.println(readName + "\tsplit-inconsistent:true");
                        }

                    } else {
                        readsToAnnotate.add(new SAMReadPair(lookup.get(readName), record));
                        ReadAnnotation ra = new ReadAnnotation();
                        ra.extractReadIntervals(lookup.get(readName), record);
                        ra.extractReadAlignmentStartEnd(lookup.get(readName), record);
                        if (readName.equals("19216203")) {
                            System.out.println();
                        }
                        if (!ra.areReadsConsistent()) {
                            System.out.println(readName + "\tsplit-inconsistent:true");
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
        pcrIndex = new HashMap<>();

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
