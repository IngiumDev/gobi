package runners;

import gtf.structs.Exon;
import gtf.structs.Interval;
import gtf.types.StrandDirection;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GTFParser;
import parsers.GenomeSequenceExtractor;
import readsimulator.ReadSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import static java.lang.System.exit;

public class ReadSimulatorRunner {
    // TODO remove refactoring, timing and logging
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("ReadSimulatorRunner").build().defaultHelp(true)
                .description("Run ReadSimulatorRunner");
        parser.addArgument("-length").required(true).help("readsimulator.Read length").metavar("<read length>").type(Integer.class);
        parser.addArgument("-frlength").required(true).help("Mean fragment length").metavar("<mean fragment length>").type(Integer.class);
        parser.addArgument("-SD").required(true).help("Standard deviation of fragment length").metavar("<standard deviation>").type(Integer.class);
        parser.addArgument("-readcounts").required(true).help("Path to read counts file").metavar("<read counts file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-mutationrate").required(true).help("Mutation rate in percent. Between 0.0 and 100.0").metavar("<mutation rate>").type(Double.class);
        parser.addArgument("-fasta").required(true).help("Path to genome FASTA file").metavar("<genome FASTA file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-fidx").required(true).help("Path to genome FASTA file index").metavar("<genome FASTA file index>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-gtf").required(true).help("Path to annotation file").metavar("<annotation file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-od").required(true).help("Output directory").metavar("<output directory>");
        if (args.length == 0) {
            parser.printHelp();
            exit(1);
        }
        try {
            Namespace res = parser.parseArgs(args);
            start(res);
        } catch (ArgumentParserException e) {
            parser.printHelp();
        }
    }

    private static void start(Namespace res) {
        long totalStartTime = System.currentTimeMillis();
        long startTime = System.currentTimeMillis();
        ReadSimulator readSimulator = new ReadSimulator.Builder()
                .setReadLength(res.getInt("length"))
                .setMeanFragmentLength(res.getInt("frlength"))
                .setFragmentLengthStandardDeviation(res.getInt("SD"))
                .setMutationRate(res.getDouble("mutationrate")/100)
                .setGtfAnnotation(GTFParser.parseGTF(res.getString("gtf")))
                .setGenomeSequenceExtractor(new GenomeSequenceExtractor(res.getString("fasta"), res.getString("fidx")))
                .setReadCounts(readCountsFile(res.getString("readcounts")))
                .build();
        
        System.out.println("Initialization time\t" + (System.currentTimeMillis() - startTime));
        startTime = System.currentTimeMillis();
        readSimulator.simulateAndWriteReads(res.getString("od"));
        System.out.println("Total Read Simulation Time\t" + (System.currentTimeMillis() - startTime));
        System.out.println("Total time\t" + (System.currentTimeMillis() - totalStartTime));
    }


    public static TreeSet<Interval> getCoveredRegion(TreeSet<Exon> exons, Interval localRegion, StrandDirection strandDirection) {
        TreeSet<Interval> coveredRegions = new TreeSet<>();
        int localStart = localRegion.getStart();
        int localEnd = localRegion.getEnd();
        int currentPos = 0;
        for (Exon exon : exons) {
            int regionLength = exon.getInterval().getLength();

            if (currentPos + regionLength > localStart) {
                int start = Math.max(exon.getInterval().getStart(), exon.getInterval().getStart() + (localStart - currentPos));
                int end = Math.min(exon.getInterval().getEnd(), exon.getInterval().getStart() + (localEnd - currentPos));
                coveredRegions.add(new Interval(start, end));
                if (currentPos + regionLength > localEnd) {
                    break;
                }
            }

            currentPos += regionLength;
        }

        return coveredRegions;
    }

    public static Map<String, Map<String, Integer>> readCountsFile(String readCountsPath) {
        Map<String, Map<String, Integer>> geneTranscriptCounts = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(readCountsPath))) {
            String line = br.readLine(); // skip first line
            while ((line = br.readLine()) != null) {
                int length = line.length();
                int tabOnePos = -1;
                int tabTwoPos = -1;
                for (int i = 0; i < length; i++) {
                    if (line.charAt(i) == '\t') {
                        if (tabOnePos == -1) {
                            tabOnePos = i;
                        } else {
                            tabTwoPos = i;
                            break;
                        }
                    }
                }
                if (tabOnePos == -1 || tabTwoPos == -1) {
                    continue; // Skip malformed lines
                }

                String geneId = line.substring(0, tabOnePos);
                String transcriptId = line.substring(tabOnePos + 1, tabTwoPos);
                int count = Integer.parseInt(line.substring(tabTwoPos + 1).trim());

                geneTranscriptCounts
                        .computeIfAbsent(geneId, k -> new HashMap<>())
                        .put(transcriptId, count);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return geneTranscriptCounts;
    }
}
