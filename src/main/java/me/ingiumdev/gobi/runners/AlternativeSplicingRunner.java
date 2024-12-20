package me.ingiumdev.gobi.runners;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.ingiumdev.gobi.bamfeatures.ReadAnnotation;
import me.ingiumdev.gobi.bamfeatures.SAMReadPair;
import me.ingiumdev.gobi.gtf.ExonSkip;
import me.ingiumdev.gobi.gtf.GTFAnnotation;
import me.ingiumdev.gobi.gtf.structs.Gene;
import me.ingiumdev.gobi.gtf.structs.Transcript;
import me.ingiumdev.gobi.gtf.treecollections.IntervalTreeForestManager;
import me.ingiumdev.gobi.gtf.treecollections.StrandUnspecificForest;
import me.ingiumdev.gobi.parsers.GTFParser;
import me.ingiumdev.gobi.readsimulator.Pair;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import static me.ingiumdev.gobi.bamfeatures.FeatureCountAnnotator.isValidRead;

public class AlternativeSplicingRunner {
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("Calculate the PSI values for exon skipping events").build().defaultHelp(true)
                .description("Run BAMFeaturesRunner");
        parser.addArgument("-gtf").required(true).help("GTF file").metavar("<gtf_file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-bam").required(true).help("BAM file").metavar("<bamfile>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-o").required(true).help("Output file").metavar("<output_tsv>");
        if (args.length == 0) {
            parser.printHelp();
            System.exit(1);
        }
        try {
            Namespace res = parser.parseArgs(args);
            start(res);
        } catch (ArgumentParserException e) {
            parser.printHelp();
        }
    }

    public static void start(Namespace res) {
        // TODO BAI
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(res.getString("bam")));
        Iterator<SAMRecord> it = reader.iterator();
        GTFAnnotation gtfAnnotation = GTFParser.parseGTF(String.valueOf(res.getString("gtf")));
        gtfAnnotation.getGenes().values().parallelStream().forEach(Gene::processIntrons);
        long startTime = System.currentTimeMillis();
        Map<String, List<ExonSkip>> exonSkips = ExonSkip.findExonSkippingEventsByGene(gtfAnnotation);
        System.out.println("LOG: Total time to find exon skipping events: " + (System.currentTimeMillis() - startTime) + " ms");

        HashMap<String, SAMRecord> lookup = new HashMap<>();
        List<SAMReadPair> readsToAnnotate = new ArrayList<>();
        String currentChromosome = "_";

        IntervalTreeForestManager forestManager = new StrandUnspecificForest();
        forestManager.init(gtfAnnotation);
        String referenceName;
        String readName;
        SAMRecord record;
        while (it.hasNext()) {
            record = it.next();
            referenceName = record.getReferenceName();
            if (!currentChromosome.equals(referenceName)) {
                lookup = new HashMap<>();
                // process reads
                List<ReadAnnotation> reads = new ArrayList<>();
                for (SAMReadPair pair : readsToAnnotate) {
                    SAMRecord first = pair.getFirst();
                    SAMRecord second = pair.getSecond();
                    ReadAnnotation readAnnotation = new ReadAnnotation(first.getReadName());
                    readAnnotation.extractReadAlignmentStartEnd(first, second);
                    List<Gene> resultGenes = forestManager.getGenesThatInclude(readAnnotation);
                    readAnnotation.extractReadIntervals(first, second);
                    readAnnotation.calculateSplit();
                    if (!resultGenes.isEmpty() && readAnnotation.areReadsConsistent()) {
                        readAnnotation.setGenesThatInclude(resultGenes);

// remove split inconsistent check if not needed
                        if (readAnnotation.findTranscriptomicMatches()) {
                            reads.add(readAnnotation);
                        }
                    }
                }
                for (ReadAnnotation read : reads) {
                    // Iterate over all genes associated with the read
                    for (Pair<Gene, List<Transcript>> match : read.getTranscriptomicMatches()) {
                        // Get the skipped exons for the current gene
                        Gene gene = match.getFirst();
                        Set<String> matchTranscriptIds = match.getSecond().stream()
                                .map(Transcript::getTranscriptID)
                                .collect(Collectors.toSet());

                    }
                }
                exonSkips.values().stream().flatMap(Collection::stream).forEach(exonSkip -> System.out.println(exonSkip.getId() + " " + exonSkip.getSV() + " " + exonSkip.getInclusionCount() + " " + exonSkip.getExclusionCount()));
                forestManager.nextTree(referenceName);
                readsToAnnotate = new ArrayList<>();
                currentChromosome = referenceName;
            }
            readName = record.getReadName();
            if (lookup.containsKey(readName)) {
                // Check if the read is the first or second of the pair
                if (record.getFirstOfPairFlag()) {
                    readsToAnnotate.add(new SAMReadPair(record, lookup.get(readName)));
                } else {
                    readsToAnnotate.add(new SAMReadPair(lookup.get(readName), record));
                }
            } else if (isValidRead(record)) {
                lookup.put(readName, record);
            }

        }
//        processLastChromosome();


    }

    private static void processNewChromosome(String referenceName) {
    }
}
