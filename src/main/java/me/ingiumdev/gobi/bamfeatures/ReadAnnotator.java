package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import me.ingiumdev.gobi.gtf.treecollections.IntervalTreeForestManager;

import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.List;

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
        return !record.getReadUnmappedFlag() && record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !ReadAnnotator.areReadsSameStrand(record) && ReadAnnotator.areReadsSameChromosome(record) && !record.getSupplementaryAlignmentFlag();
    }

    public static boolean areReadsSameChromosome(SAMRecord record) {
        return record.getReferenceName().equals(record.getMateReferenceName());
    }

    public abstract void annotateReads();

    protected abstract void processNewChromosome(String referenceName);

    protected abstract ReadAnnotation processRead(SAMReadPair samReadPair);

    protected abstract void processLastChromosome();
}
