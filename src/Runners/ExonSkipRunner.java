package Runners;

import gtf.*;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.helper.HelpScreenException;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GTFParser;

import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;


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

        // try to find exon with multiple cds
/* given a SV-intron candidate we have
to check if there exists at least one WT
Given Intron I as SV candidate:
 The wild types have two introns:
1) one starting at I.getX1() -
2) one ending at I.getX2()
3) but not at the same time
4)(that would be the SV-intron)
→ if we calculate the set of isoforms having such introns for all three cases,
we can simply derive the set of corresponding WT isoforms (or that this is empty)
Calculate the following sets Set<String> (entries: cds ids):

the set SV (CDS-s containing I)

the set WT_start (CDS-s having introns starting I.getX1())

the set WT_end (CDS-s having introns ending I.getX2())
 the set WT= (WT_start intersect WT_end) – SV
→ if WT is not empty we have an exon skipping event.
• id (gene id)
• symbol (gene symbol)
• chr (chromosome)
• strand (+ or -)
• nprots (number of annotated CDS in the gene)
• ntrans (number of annotated transcripts in the gene)
• SV (the SV intron as start:end)
• WT (all the WT introns within the SV intron separated by | as start:end)
• SV prots (ids of the SV CDS, separated by |)
• WT prots (ids of the WT CDS, separated by |)
• min skipped exon the minimal number of skipped exons in any WT/SV pair
• max skipped exon the maximum number of skipped exons in any WT/SV pair
• min skipped bases the minimal number of skipped bases (joint length of skipped exons)
in any WT/SV pair
• max skipped bases the maximum number of skipped bases (joint length of skipped
exons) in any WT/SV pair.

*/
        // Iterate through each gene, then through each intron
        Set<ExonSkip> exonSkips = new HashSet<>();
        // Iterate through each gene, then for each gene iterate through each transcript so that we can iterate through the introns
        for (Gene gene : annotation.getGenes().values()) {
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
                            String protein_id = t.getCds().getFirst().getAttribute("protein_id");
                            if (i.equals(intron)) {
                                SV.add(protein_id);
                            } else if (i.getStart() == intron.getStart()) {
                                WT_start.add(protein_id);
                                WT_introns.add(i);
                            } else if (i.getEnd() == intron.getEnd()) {
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
                        exonSkip.setSymbol(gene.getAttribute("gene_name"));
                        exonSkip.setChr(gene.getSeqname());
                        exonSkip.setStrand(gene.getStrand());
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
        }

        System.out.println("done");
    }


}
