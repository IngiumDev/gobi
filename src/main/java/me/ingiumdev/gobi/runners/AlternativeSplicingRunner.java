package me.ingiumdev.gobi.runners;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.ingiumdev.gobi.bamfeatures.AlternativeSplicingAnnotator;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.File;

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
        SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .setOption(SamReaderFactory.Option.EAGERLY_DECODE, false)
                .open(new File(res.getString("bam")));

        AlternativeSplicingAnnotator alternativeSplicingAnnotator = new AlternativeSplicingAnnotator.Builder()
                .setSamReader(reader)
                .setGtfFile(new File(res.getString("gtf")))
                .setOutputFile(new File(res.getString("o")))
                .build();
        alternativeSplicingAnnotator.init();
        alternativeSplicingAnnotator.annotateReads();
    }

    private static void processNewChromosome(String referenceName) {
    }
}
