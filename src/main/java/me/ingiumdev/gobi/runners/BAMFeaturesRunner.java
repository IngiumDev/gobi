package me.ingiumdev.gobi.runners;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.ingiumdev.gobi.bamfeatures.ReadAnnotator;
import me.ingiumdev.gobi.gtf.types.StrandDirection;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.File;

public class BAMFeaturesRunner {

    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("BAMFeaturesRunner").build().defaultHelp(true)
                .description("Run BAMFeaturesRunner");
        parser.addArgument("-gtf").required(true).help("GTF file").metavar("<gtf_file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-bam").required(true).help("BAM file").metavar("<bamfile>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-o").required(true).help("Output file").metavar("<output_tsv>");
        parser.addArgument("-frstrand").help("true/false").metavar("<true/false>");
        parser.addArgument("-analysis").help("Path to the analysis file").metavar("<analysis-file-path>");
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

    // TODO: Ignore unneeded chromosomes, ignore cds, etc
    private static void start(Namespace res) {
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(res.getString("bam")));
        StrandDirection strandSpecific = (res.getString("frstrand") == null) ? StrandDirection.UNSPECIFIED : (res.getString("frstrand").equals("true") ? StrandDirection.FORWARD : StrandDirection.REVERSE);
        ReadAnnotator readAnnotator = new ReadAnnotator.Builder()
                .setSamReader(reader)
                .setGtfFile(new File(res.getString("gtf")))
                .setOutputFile(new File(res.getString("o")))
                .setStrandSpecificity(strandSpecific)
                .build();
        long start = System.currentTimeMillis();
        readAnnotator.annotateReads();

        System.out.println("Time annotate: " + (System.currentTimeMillis() - start) + "ms");
        // todo close samreader
    }
}
