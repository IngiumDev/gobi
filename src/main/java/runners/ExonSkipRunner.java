package runners;

import gtf.ExonSkip;
import gtf.GTFAnnotation;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GTFParser;

import java.util.List;


public class ExonSkipRunner {
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("ExonSkipRunner").build().defaultHelp(true)
                .description("Run ExonSkipRunner");
        parser.addArgument("-gtf").required(true).help("GTF file").metavar("<GTF file>");
        parser.addArgument("-o").required(true).help("Output file").metavar("<output file path>");

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
        GTFAnnotation GTFAnnotation = GTFParser.parseGTF(res.getString("gtf"));
        long start = System.currentTimeMillis();
        List<ExonSkip> exonSkips = ExonSkip.findExonSkippingEvents(GTFAnnotation);
        long end = System.currentTimeMillis();
        System.out.println("Time taken to process exon skipping events: " + (end - start) + "ms");
        // Write the exon skipping events to a file
        ExonSkip.writeExonSkipToFile(res.getString("o"), exonSkips);


    }


}
