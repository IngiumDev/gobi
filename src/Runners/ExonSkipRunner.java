package Runners;

import gtf.*;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.helper.HelpScreenException;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GTFParser;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
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
        System.out.println("Parsing GTF file");
        Annotation annotation = GTFParser.parseGTF(res.getString("gtf"));
        // Calculate the set of introns for each gene
        System.out.println("Calculating introns");

        // Iterate through each gene, then through each intron

        Set<ExonSkip> exonSkips = findExonSkippingEvens(annotation);

        // Write the exon skipping events to a file
        System.out.println("Writing exon skipping events to file");
        writeExonSkipToFile(res.getString("o"), exonSkips);


        System.out.println("done");
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
                                exonSkip.getWT_prots().stream().collect(Collectors.joining("|")) + "\t" +
                                exonSkip.getSV_prots().stream().collect(Collectors.joining("|")) + "\t" +
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

    private static Set<ExonSkip> findExonSkippingEvens(Annotation annotation) {
        Set<ExonSkip> exonSkips = new HashSet<>();
        annotation.getGenes().values().parallelStream().forEach(gene -> {
            for (Transcript transcript : gene.getTranscripts().values()) {
                for (Interval intron : transcript.getIntrons()) {
                    // Check if the intron is a candidate for exon skipping
                    Set<String> SV = new HashSet<>();
                    Set<String> WT_start = new HashSet<>();
                    Set<String> WT_end = new HashSet<>();
                    Set<String> WT = new HashSet<>();
                    Set<Interval> WT_introns = new HashSet<>();

                    int minSkippedExons = Integer.MAX_VALUE;
                    int maxSkippedExons = Integer.MIN_VALUE;
                    int minSkippedBases = Integer.MAX_VALUE;
                    int maxSkippedBases = Integer.MIN_VALUE;

                    for (Transcript t : gene.getTranscripts().values()) {
                        // only consider transcripts with cds
                        if (t.getCds().size() == 0) {
                            //System.out.println("No CDS found for transcript " + t.getId());
                            continue;
                        }

                        for (Interval i : t.getIntrons()) {
                            // Check if the intron is the same as the candidate intron (SV)
                            if (i.equals(intron)) {
                                String protein_id = t.getCds().getFirst().getAttribute("protein_id");
                                SV.add(protein_id);
                            } else if (i.getStart() == intron.getStart()) {
                                String protein_id = t.getCds().getFirst().getAttribute("protein_id");
                                WT_start.add(protein_id);
                                WT_introns.add(i);
                            } else if (i.getEnd() == intron.getEnd()) {
                                String protein_id = t.getCds().getFirst().getAttribute("protein_id");
                                WT_end.add(protein_id);
                                WT_introns.add(i);
                            }
                        }
                    }
                    WT.addAll(WT_start);
                    WT.retainAll(WT_end);
                    WT.removeAll(SV);
                    if (!WT.isEmpty()) {
                        // We have an exon skipping event
                        // TODO add exon skip to the list
                        ExonSkip exonSkip = new ExonSkip();
                        exonSkip.setId(gene.getId());
                        // If gene is directly annoted
                        // TODO: Migrate to moving data up the chain from exon to transcript to gene
                        if (gene.getAttributes() != null) {

                            exonSkip.setSymbol(gene.getAttribute("gene_name"));
                            exonSkip.setChr(gene.getSeqname());
                            exonSkip.setStrand(gene.getStrand());
                        } else {
                            // TODO Migrate away from storing strand in anything other than gene
                            Exon firstExon = transcript.getExons().first();
                            exonSkip.setSymbol(firstExon.getAttribute("gene_name"));
                            exonSkip.setChr(firstExon.getSeqname());
                            exonSkip.setStrand(firstExon.getStrand());
                        }

                        // nprots is the number of transcripts where the cds is not empty
                        exonSkip.setNprots((int) gene.getTranscripts().values().stream().filter(t -> !t.getCds().isEmpty()).count());
                        exonSkip.setNtrans(gene.getTranscripts().size());
                        exonSkip.setSV(new Interval(intron));
                        exonSkip.setWT(new TreeSet<>(WT_introns));
                        exonSkip.setSV_prots(SV);
                        exonSkip.setWT_prots(WT);
                        // Calculate min and max skipped exons and bases

                        for (String protein_id : WT) {
                            // Get the transcript
                            String transcript_id = gene.convertProteinToTranscript(protein_id);
                            Transcript t = gene.getTranscript(transcript_id);
                            int skippedExons = 0;
                            int skippedBases = 0;
                            // Find the number of exons that lie within the intron
                            for (Exon exon : t.getExons()) {
                                if (intron.getStart() <= exon.getInterval().getStart() && intron.getEnd() >= exon.getInterval().getEnd()) {
                                    skippedExons++;
                                    skippedBases += exon.getInterval().getLength();
                                }
                            }
                            minSkippedExons = Math.min(minSkippedExons, skippedExons);
                            maxSkippedExons = Math.max(maxSkippedExons, skippedExons);
                            minSkippedBases = Math.min(minSkippedBases, skippedBases);
                            maxSkippedBases = Math.max(maxSkippedBases, skippedBases);

                        }

                        exonSkip.setMax_skipped_exon(maxSkippedExons);
                        exonSkip.setMin_skipped_exon(minSkippedExons);
                        exonSkip.setMax_skipped_bases(maxSkippedBases);
                        exonSkip.setMin_skipped_bases(minSkippedBases);

                        exonSkips.add(exonSkip);
                    }
                }
            }
        });
        return exonSkips;
    }


}
