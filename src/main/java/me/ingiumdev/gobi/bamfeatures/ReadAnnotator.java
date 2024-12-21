package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import me.ingiumdev.gobi.gtf.treecollections.IntervalTreeForestManager;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public abstract class ReadAnnotator {
    protected SamReader samReader;
    protected File outputFile;
    protected File gtfFile;
    protected List<SAMReadPair> readsToAnnotate;
    protected BufferedWriter writer;
    protected IntervalTreeForestManager forestManager;
    protected HashMap<String, SAMRecord> lookup;
    protected String currentChromosome = "_";

    public ReadAnnotator() {
    }

    public static boolean areReadsSameStrand(SAMRecord record) {
        return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
    }

    public static boolean isValidRead(SAMRecord record) {
        return !record.getReadUnmappedFlag() && record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !ReadAnnotator.areReadsSameStrand(record) && !record.getNotPrimaryAlignmentFlag() && ReadAnnotator.areReadsSameChromosome(record) && !record.getSupplementaryAlignmentFlag();
    }

    public static boolean areReadsSameChromosome(SAMRecord record) {
        return record.getReferenceName().equals(record.getMateReferenceName());
    }

    public void annotateReads() {
        long startTime = System.currentTimeMillis();
        Iterator<SAMRecord> it = samReader.iterator();
        String referenceName;
        String readName;
        SAMRecord record;
        // TODO: remove from lookup once it's processed
        Set<String> seennames = new HashSet<>();
        while (it.hasNext()) {
            record = it.next();
            referenceName = record.getReferenceName();
            if (!currentChromosome.equals(referenceName)) {
                processNewChromosome(referenceName);
            }
            readName = record.getReadName();
            if (lookup.containsKey(readName) && isValidRead(record)) {
                // Check if the read is the first or second of the pair
                boolean seen = seennames.add(readName);
                if (true) {
                    if (record.getFirstOfPairFlag()) {
                        readsToAnnotate.add(new SAMReadPair(record, lookup.get(readName)));
                    } else {
                        readsToAnnotate.add(new SAMReadPair(lookup.get(readName), record));
                    }
                } else {
                    if (readName.equals("77779")) {
                        System.out.println("Duplicate read found");
                    }
                }
            } else if (isValidRead(record)) {
                lookup.put(readName, record);
            }
        }
        processLastChromosome();
        System.out.println("Annotation Time taken: " + (System.currentTimeMillis() - startTime) + "ms");
    }

    public abstract void init();

    protected abstract void processNewChromosome(String referenceName);

    protected abstract ReadAnnotation processRead(SAMReadPair samReadPair);

    protected abstract void processLastChromosome();
}
