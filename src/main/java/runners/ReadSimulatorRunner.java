package runners;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ReadSimulatorRunner {
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("ReadSimulatorRunner").build().defaultHelp(true)
                .description("Run ReadSimulatorRunner");

        parser.addArgument("-length").required(true).help("Read length").metavar("<read length>").type(Integer.class);
        parser.addArgument("-frlength").required(true).help("Mean fragment length").metavar("<mean fragment length>").type(Integer.class);
        parser.addArgument("-SD").required(true).help("Standard deviation of fragment length").metavar("<standard deviation>").type(Integer.class);
        parser.addArgument("-readcounts").required(true).help("Path to read counts file").metavar("<read counts file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-mutationrate").required(true).help("Mutation rate in percent. Between 0.0 and 100.0").metavar("<mutation rate>").type(Double.class);
        parser.addArgument("-fasta").required(true).help("Path to genome FASTA file").metavar("<genome FASTA file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-fidx").required(true).help("Path to genome FASTA file index").metavar("<genome FASTA file index>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-gtf").required(true).help("Path to annotation file").metavar("<annotation file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-od").required(true).help("Output directory").metavar("<output directory>").type(new FileArgumentType().verifyIsDirectory());
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

    private static void start(Namespace res) {
    }
}
