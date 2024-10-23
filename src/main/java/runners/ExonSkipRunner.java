package runners;

import gtf.*;
import gtf.structs.Interval;
import gtf.types.StrandDirection;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.helper.HelpScreenException;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GTFParser;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.stream.Collectors;


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
        GTFAnnotation GTFAnnotation = GTFParser.parseGTF(res.getString("gtf"));

        // Iterate through each gene, then through each intron
        long start = System.currentTimeMillis();
        Set<ExonSkip> exonSkips = ExonSkip.findExonSkippingEvents(GTFAnnotation);
        long end = System.currentTimeMillis();
        System.out.println("Time taken: " + (end - start) + "ms");
        // Write the exon skipping events to a file
        writeExonSkipToFile(res.getString("o"), exonSkips);


    }

    public static void writeExonSkipToFile(String o, Set<ExonSkip> exonSkips) {

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(o))) {
            writer.write("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tSV_prots\tWT_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n");
            for (ExonSkip exonSkip : exonSkips) {
                writer.write(
                        exonSkip.getId() + "\t" +
                                exonSkip.getSymbol() + "\t" +
                                exonSkip.getChr() + "\t" +
                                (exonSkip.getStrand() == StrandDirection.FORWARD ? "+" : "-") + "\t" +
                                exonSkip.getNprots() + "\t" +
                                exonSkip.getNtrans() + "\t" +
                                exonSkip.getSV() + "\t" +
                                exonSkip.getWT().stream().map(Interval::toString).collect(Collectors.joining("|")) + "\t" +
                                String.join("|", exonSkip.getWT_prots()) + "\t" +
                                String.join("|", exonSkip.getSV_prots()) + "\t" +
                                exonSkip.getMin_skipped_exon() + "\t" +
                                exonSkip.getMax_skipped_exon() + "\t" +
                                exonSkip.getMin_skipped_bases() + "\t" +
                                exonSkip.getMax_skipped_bases() + "\n"
                );
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }


}
