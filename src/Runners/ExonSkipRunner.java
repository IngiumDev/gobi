package Runners;

import gtf.Annotation;
import parsers.GTFParser;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.helper.HelpScreenException;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.Namespace;

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
        } catch (HelpScreenException e) {
        parser.printHelp();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void start(Namespace res) {
        Annotation annotation = GTFParser.parseGTF(res.getString("gtf"));
        System.out.println();
    }
}
