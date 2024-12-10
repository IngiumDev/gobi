package runners;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;


public class GeneLengthRunner {
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("GeneLengthRunner").build().defaultHelp(true)
                .description("Run GeneLengthRunner");
        parser.addArgument("-gtf").required(true).help("GTF file").metavar("<gtf_file>").type(new FileArgumentType().verifyIsFile());
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
        // We want to parse a gtf and output a tsv with gene id gene name chromosome start end length

    }
}
