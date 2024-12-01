package runners;

import gtf.structs.Exon;
import gtf.structs.Interval;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GenomeSequenceExtractor;
import readsimulator.ReadSimulator;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;

import static java.lang.System.exit;

public class ReadSimulatorRunner {
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
        parser.addArgument("-analysis-script").help("Path to the analysis file to analyze the results directly").metavar("<path-to-script>").type(new FileArgumentType().verifyIsFile());
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
                .setReadCounts(readCountsFile(res.getString("readcounts")))
                .setGtfAnnotation(res.getString("gtf"))
                .setGenomeSequenceExtractor(new GenomeSequenceExtractor(res.getString("fasta"), res.getString("fidx")))
                .build();
        System.out.println("Initialization time\t" + (System.currentTimeMillis() - startTime));
        startTime = System.currentTimeMillis();
        readSimulator.simulateAndWriteReads(res.getString("od"));
        System.out.println("Total Read Simulation Time\t" + (System.currentTimeMillis() - startTime));
        System.out.println("Total time\t" + (System.currentTimeMillis() - totalStartTime));
        if (res.get("analysis_script") != null) {
            String scriptPath = ((File) res.get("analysis_script")).getAbsolutePath();

            // Use ProcessBuilder to call the Python script
            List<String> command = new ArrayList<>();
            command.add("python3"); // or "python3", depending on your environment
            command.add(scriptPath);
            command.add(Paths.get(res.getString("od"), "read.mappinginfo").toString());
            command.add(res.getString("od"));

            ProcessBuilder pb = new ProcessBuilder(command);
            pb.redirectErrorStream();

            try {
                Process process = pb.start();
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line); // Print the output
                }

                // Wait for the process to complete and get the exit code
                int exitCode = process.waitFor();
                if (exitCode == 0) {
                    System.out.println("Python script executed successfully.");
                } else {
                    System.err.println("Python script execution failed with exit code: " + exitCode);
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }

        }


    }


    public static TreeSet<Interval> getCoveredRegion(TreeSet<Exon> exons, Interval localRegion) {
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
    public static TreeSet<Interval> cut(TreeSet<Exon> exons, Interval cutRegion) {
        TreeSet<Interval> cutRegions = new TreeSet<>();
        int cutStart = cutRegion.getStart();
        int cutEnd = cutRegion.getEnd();

        for (Exon exon : exons) {
            Interval exonInterval = exon.getInterval();

            // Check if the exon overlaps with the cut region
            if (exonInterval.getEnd() >= cutStart && exonInterval.getStart() <= cutEnd) {
                // Calculate the overlap
                int start = Math.max(exonInterval.getStart(), cutStart);
                int end = Math.min(exonInterval.getEnd(), cutEnd);

                // Add the overlapping region
                cutRegions.add(new Interval(start, end));
            }

            // If the current exon ends beyond the cut region, stop processing
            if (exonInterval.getStart() > cutEnd) {
                break;
            }
        }

        return cutRegions;
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
