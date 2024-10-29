package runners;

import gtf.ExonSkip;
import gtf.GTFAnnotation;
import gtf.structs.GTFTimer;
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
        parser.addArgument("-a", "--analysis").required(false).help("File Path to the analysis file").metavar("<analysis file path>");
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
        // GTF Parsing
        long totalStartTime = System.currentTimeMillis();
        GTFAnnotation GTFAnnotation = GTFParser.parseGTF(res.getString("gtf"));

        // Exon Skipping Calculation
        long startTime = System.currentTimeMillis();
        List<ExonSkip> exonSkips = ExonSkip.findExonSkippingEvents(GTFAnnotation);
        GTFTimer.setExonProcessTime((System.currentTimeMillis() - startTime));
        System.out.println("LOG: Total time to find exon skipping events: " + GTFTimer.getExonProcessTime() + " ms");

        // Output to file
        startTime = System.currentTimeMillis();
        ExonSkip.writeExonSkipToFile(res.getString("o"), exonSkips);
        GTFTimer.setOutputTime(System.currentTimeMillis() - startTime);
        System.out.println("LOG: Total time to write exon skipping events: " + GTFTimer.getOutputTime() + " ms");
        GTFTimer.setTotalTime(System.currentTimeMillis() - totalStartTime);
        System.out.println("LOG: Total time: " + GTFTimer.getTotalTime() + " ms");

        // Analysis (optional)
        if (res.getString("analysis") != null) {
            startTime = System.currentTimeMillis();
            ExonSkip.analyzeExonSkippingEvents(GTFAnnotation, exonSkips, res.getString("analysis"), res.getString("gtf"));
            System.out.println("LOG: Total time to analyze exon skipping events: " + (System.currentTimeMillis() - startTime) + " ms");
        }
    }


}
