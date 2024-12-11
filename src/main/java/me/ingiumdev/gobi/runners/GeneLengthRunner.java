package me.ingiumdev.gobi.runners;

import me.ingiumdev.gobi.gtf.GTFAnnotation;
import me.ingiumdev.gobi.gtf.structs.Exon;
import me.ingiumdev.gobi.gtf.structs.Gene;
import me.ingiumdev.gobi.gtf.structs.Interval;
import me.ingiumdev.gobi.gtf.structs.Transcript;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.TreeSet;

import static me.ingiumdev.gobi.bamfeatures.ReadAnnotation.mergeInterval2;
import static me.ingiumdev.gobi.parsers.GTFParser.parseGTF;


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
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(res.getString("o")))) {
            bw.write("gene_id\tgene_name\tchromosome\tstart\tend\tlength");
            GTFAnnotation gtf = parseGTF(res.getString("gtf"));
            for (String geneId : gtf.getGenes().keySet()) {
                Gene gene = gtf.getGenes().get(geneId);
                int mergedLength = calculateMergedGeneLength(gene);
                bw.newLine();
                bw.write(geneId + "\t" + gene.getGeneName() + "\t" + gene.getSeqname() + "\t" + gene.getStart() + "\t" + gene.getStop() + "\t" + mergedLength);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static int calculateMergedGeneLength(Gene gene) {
        TreeSet<Interval> mergedGeneTranscript = new TreeSet<>();
        for (Transcript transcript : gene.getTranscripts().values()) {
            for (Exon exon : transcript.getExons()) {
                mergeInterval2(mergedGeneTranscript, exon.getInterval());
            }
        }
        return mergedGeneTranscript.stream()
                .mapToInt(interval -> interval.getEnd() - interval.getStart() + 1)
                .sum();
    }
}
