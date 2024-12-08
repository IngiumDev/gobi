package gtf;

import gtf.structs.*;
import gtf.types.StrandDirection;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class ExonSkip {

    //gene_id
    private final String id;
    // gene symbol
    private final String symbol;
    // chromosome
    private final String chr;
    // strand
    private final StrandDirection strand;
    // number of annotated CDS in the gene
    private final int nprots;
    // number of annotated transcripts in the gene
    private final int ntrans;
    // the SV intron as start:end
    private final Interval SV;
    // all the WT introns within the SV intron separated by | as start:end
    private final TreeSet<Interval> WT;
    // ids of the SV CDS, separated by | (protein_id)
    private final Set<String> SV_prots;
    // ids of the WT CDS, separated by | (protein_id)
    private final Set<String> WT_prots;
    // the minimal number of skipped exons in any WT/SV pair
    private final int min_skipped_exon;
    // the maximum number of skipped exons in any WT/SV pair
    private final int max_skipped_exon;
    // the minimal number of skipped bases (joint length of skipped exons) in any WT/SV pair
    private final int min_skipped_bases;
    // the maximum number of skipped bases (joint length of skipped exons) in any WT/SV pair
    private final int max_skipped_bases;

    private ExonSkip(Builder builder) {
        this.id = builder.id;
        this.symbol = builder.symbol;
        this.chr = builder.chr;
        this.strand = builder.strand;
        this.nprots = builder.nprots;
        this.ntrans = builder.ntrans;
        this.SV = builder.SV;
        this.WT = builder.WT;
        this.SV_prots = builder.SV_prots;
        this.WT_prots = builder.WT_prots;
        this.min_skipped_exon = builder.min_skipped_exon;
        this.max_skipped_exon = builder.max_skipped_exon;
        this.min_skipped_bases = builder.min_skipped_bases;
        this.max_skipped_bases = builder.max_skipped_bases;
    }

    public static List<ExonSkip> findExonSkippingEvents(GTFAnnotation gtfAnnotation) {
        // TODO use .collect() to avoid synchronization
        List<ExonSkip> exonSkips = Collections.synchronizedList(new ArrayList<>());
        gtfAnnotation.getGenes().values().parallelStream().forEach(gene -> {
            for (Interval intronCandidate : gene.getIntrons()) {
                Set<String> spliceVariantTranscripts = new HashSet<>();
                Set<String> wildTypeTranscripts = new HashSet<>();
                Set<Interval> wildTypeIntrons = new TreeSet<>();
                checkIntronCandidate(gene, intronCandidate, spliceVariantTranscripts, wildTypeTranscripts, wildTypeIntrons);

                if (!spliceVariantTranscripts.isEmpty() && !wildTypeTranscripts.isEmpty()) {
                    ExonSkip exonSkip = createExonSkipEvent(gene, intronCandidate, spliceVariantTranscripts, wildTypeTranscripts, wildTypeIntrons);
                    exonSkips.add(exonSkip);
                }
            }
        });
        return exonSkips;
    }

    private static void checkIntronCandidate(Gene gene, Interval intronCandidate, Set<String> spliceVariantTranscripts, Set<String> wildTypeTranscripts, Set<Interval> wildTypeIntrons) {
        for (Transcript transcriptToCheck : gene.getTranscripts().values()) {
            Set<Interval> intronsToAdd = new TreeSet<>();

            int intronCandidateStart = intronCandidate.getStart();
            int intronCandidateEnd = intronCandidate.getEnd();
            boolean isStartIntron = false;
            boolean isEndIntron = false;
            CodingSequence previousCds = null;


            for (CodingSequence cds : transcriptToCheck.getCds()) {
                if (previousCds != null) {
                    if (previousCds.getInterval().getEnd() + 1 >= intronCandidateStart && cds.getInterval().getStart() - 1 <= intronCandidateEnd) {
                        if (previousCds.getInterval().getEnd() + 1 == intronCandidateStart && cds.getInterval().getStart() - 1 == intronCandidateEnd) {
                            spliceVariantTranscripts.add(transcriptToCheck.getTranscriptID());
                            break;
                        } else if (previousCds.getInterval().getEnd() + 1 == intronCandidateStart) {
                            isStartIntron = true;
                            intronsToAdd.add(new Interval(intronCandidateStart, cds.getInterval().getStart() - 1));
                        } else if (cds.getInterval().getStart() - 1 == intronCandidateEnd) {
                            isEndIntron = true;
                            intronsToAdd.add(new Interval(previousCds.getInterval().getEnd() + 1, intronCandidateEnd));
                        } else {
                            intronsToAdd.add(new Interval(previousCds.getInterval().getEnd() + 1, cds.getInterval().getStart() - 1));
                        }
                    }
                }
                previousCds = cds;
            }
            if (isStartIntron && isEndIntron) {
                wildTypeTranscripts.add(transcriptToCheck.getTranscriptID());
                wildTypeIntrons.addAll(intronsToAdd);
            }
        }

    }

    private static ExonSkip createExonSkipEvent(Gene gene, Interval intronCandidate, Set<String> spliceVariantTranscripts, Set<String> wildTypeTranscripts, Set<Interval> wildTypeIntrons) {
        SkippedExonsBases skippedExonsBases = calculateSkippedExonsAndBases(gene, intronCandidate, wildTypeTranscripts);
        Set<String> spliceVariantProteins = convertTranscriptToProteinID(gene, spliceVariantTranscripts);
        Set<String> wildTypeProteins = convertTranscriptToProteinID(gene, wildTypeTranscripts);
        return new Builder()
                .setId(gene.getGeneID())
                .setSymbol(gene.getGeneName())
                .setChr(gene.getSeqname())
                .setStrand(gene.getStrand())
                .setNprots((int) gene.getTranscripts().values().stream().map(Transcript::getCds).filter(cds -> cds != null && !cds.isEmpty()).count())
                .setNtrans(gene.getTranscripts().size())
                .setSV(new Interval(intronCandidate))
                .setWT(new TreeSet<>(wildTypeIntrons))
                .setSV_prots(spliceVariantProteins)
                .setWT_prots(wildTypeProteins)
                .setMax_skipped_exon(skippedExonsBases.maxSkippedExons())
                .setMin_skipped_exon(skippedExonsBases.minSkippedExons())
                .setMax_skipped_bases(skippedExonsBases.maxSkippedBases())
                .setMin_skipped_bases(skippedExonsBases.minSkippedBases())
                .build();
    }

    private static Set<String> convertTranscriptToProteinID(Gene gene, Set<String> transcriptIDs) {
        return transcriptIDs.stream()
                .map(gene::getTranscript)
                .map(transcript -> {
                    String proteinID = transcript.getCds().getFirst().getProteinID();
                    if (proteinID != null) {
                        return proteinID;
                    }
                    String ccdsID = transcript.getCds().getFirst().getCcdsID();
                    return Objects.requireNonNullElse(ccdsID, "NONE");
                }).collect(Collectors.toSet());
    }

    private static SkippedExonsBases calculateSkippedExonsAndBases(Gene gene, Interval intron, Set<String> wildTypeTranscripts) {
        int minSkippedExons = Integer.MAX_VALUE;
        int maxSkippedExons = Integer.MIN_VALUE;
        int minSkippedBases = Integer.MAX_VALUE;
        int maxSkippedBases = Integer.MIN_VALUE;
        for (String transcript_id : wildTypeTranscripts) {
            // Get the transcript
            Transcript t = gene.getTranscript(transcript_id);
            int skippedExons = 0;
            int skippedBases = 0;
            // Find the number of exons that lie within the intron
            for (CodingSequence cds : t.getCds()) {
                if (intron.getStart() <= cds.getInterval().getStart() && intron.getEnd() >= cds.getInterval().getEnd()) {
                    skippedExons++;
                    skippedBases += cds.getInterval().getLength();
                }
            }
            minSkippedExons = Math.min(minSkippedExons, skippedExons);
            maxSkippedExons = Math.max(maxSkippedExons, skippedExons);
            minSkippedBases = Math.min(minSkippedBases, skippedBases);
            maxSkippedBases = Math.max(maxSkippedBases, skippedBases);

        }
        return new SkippedExonsBases(minSkippedExons, maxSkippedExons, minSkippedBases, maxSkippedBases);
    }

    public static void writeExonSkipToFile(String o, List<ExonSkip> exonSkips) {

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
                                String.join("|", exonSkip.getSV_prots()) + "\t" +
                                String.join("|", exonSkip.getWT_prots()) + "\t" +
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

    public static void analyzeExonSkippingEvents(GTFAnnotation gtfAnnotation, List<ExonSkip> exonSkips, String analysis, String gtf) {
        // What we want: first of all general statistics like how long the file is, how long the total time was, total amount of genes, total gene length, total number of transcipts, total number of transcripts with cds, total number of exons, total number of cds, total number of introns, introns per gene, cds per gene
        // Then we want to list the genes: id and name, length, number exon skipping events, number of transcripts, number of cds, number of exons, number of introns
        // General statistics
        long fileLength = -1;
        try {
            fileLength = Files.lines(Paths.get(gtf)).count();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        StringBuilder output = new StringBuilder("# General statistics\n" +
                "File length\t" + fileLength + "\n" +
                "GTF parse time\t" + GTFTimer.getGtfParseTime() + "\n" +
                "Intron process time\t" + GTFTimer.getIntronProcessTime() + "\n" +
                "Exon process time\t" + GTFTimer.getExonProcessTime() + "\n" +
                "Output time\t" + GTFTimer.getOutputTime() + "\n" +
                "Total time taken\t" + GTFTimer.getTotalTime() + "\n");
        int totalGenes = gtfAnnotation.getGenes().size();
        output.append("Total genes\t").append(totalGenes).append("\n");
        int totalGenesWithCDS = (int) gtfAnnotation.getGenes().values().stream()
                .filter(gene -> gene.getTranscripts().values().stream().anyMatch(transcript -> !transcript.getCds().isEmpty()))
                .count();
        output.append("Total genes with CDS\t").append(totalGenesWithCDS).append("\n");
        int totalGenesWithoutCDS = totalGenes - totalGenesWithCDS;

        output.append("Total genes without CDS\t").append(totalGenesWithoutCDS).append("\n");
        long totalGeneLength = gtfAnnotation.getGenes().values().stream()
                .mapToLong(gene -> {
                    Interval interval = gene.getInterval();
                    if (interval != null) {
                        return interval.getLength();
                    }
                    TreeSet<Exon> allExons = gene.getTranscripts().values().stream()
                            .flatMap(transcript -> transcript.getExons().stream())
                            .collect(Collectors.toCollection(TreeSet::new));

                    if (allExons.isEmpty()) {
                        return 0L;
                    }

                    long start = allExons.first().getInterval().getStart();
                    long end = allExons.last().getInterval().getEnd();
                    return end - start + 1;
                })
                .sum();

        output.append("Total gene length\t").append(totalGeneLength).append("\n");

        long totalGeneLengthWithCDS = gtfAnnotation.getGenes().values().stream()
                .filter(gene -> gene.getTranscripts().values().stream().anyMatch(transcript -> !transcript.getCds().isEmpty()))
                .mapToLong(gene -> {
                    Interval interval = gene.getInterval();
                    if (interval != null) {
                        return interval.getLength();
                    }
                    TreeSet<Exon> allExons = gene.getTranscripts().values().stream()
                            .flatMap(transcript -> transcript.getExons().stream())
                            .collect(Collectors.toCollection(TreeSet::new));

                    if (allExons.isEmpty()) {
                        return 0L;
                    }

                    long start = allExons.first().getInterval().getStart();
                    long end = allExons.last().getInterval().getEnd();
                    return end - start + 1;
                })
                .sum();

        output.append("Total gene length with CDS\t").append(totalGeneLengthWithCDS).append("\n");
        int totalTranscripts = gtfAnnotation.getGenes().values().stream()
                .mapToInt(gene -> gene.getTranscripts().size())
                .sum();
        output.append("Total transcripts\t").append(totalTranscripts).append("\n");


        int totalTranscriptsWithCDS = gtfAnnotation.getGenes().values().stream()
                .mapToInt(gene -> (int) gene.getTranscripts().values().stream().filter(transcript -> !transcript.getCds().isEmpty()).count())
                .sum();
        output.append("Total transcripts with CDS\t").append(totalTranscriptsWithCDS).append("\n");
        int totalTranscriptsWithoutCDs = totalTranscripts - totalTranscriptsWithCDS;
        output.append("Total transcripts without CDS\t").append(totalTranscriptsWithoutCDs).append("\n");
        int totalExons = gtfAnnotation.getGenes().values().stream()
                .flatMap(gene -> gene.getTranscripts().values().stream())
                .mapToInt(transcript -> transcript.getExons().size())
                .sum();
        output.append("Total exons\t").append(totalExons).append("\n");
        int totalCDS = gtfAnnotation.getGenes().values().stream()
                .flatMap(gene -> gene.getTranscripts().values().stream())
                .mapToInt(transcript -> transcript.getCds().size())
                .sum();
        output.append("Total CDS\t").append(totalCDS).append("\n");
        int totalIntrons = gtfAnnotation.getGenes().values().stream()
                .mapToInt(gene -> gene.getIntrons().size())
                .sum();
        output.append("Total introns\t").append(totalIntrons).append("\n");
        int totalExonSkippingEvents = exonSkips.size();
        output.append("Total exon skipping events\t").append(totalExonSkippingEvents).append("\n");
        double transcriptsPerGene = (double) totalTranscripts / totalGenes;
        output.append("Transcripts per gene\t").append(transcriptsPerGene).append("\n");
        double cdsPerGene = (double) totalCDS / totalGenes;
        output.append("CDS per gene\t").append(cdsPerGene).append("\n");
        double intronsPerGene = (double) totalIntrons / totalGenes;
        output.append("Introns per gene\t").append(intronsPerGene).append("\n");
        double exonsPerGene = (double) totalExons / totalGenes;
        output.append("Exons per gene\t").append(exonsPerGene).append("\n");


        // Gene statistics
        Map<String, Long> geneIdToTotalExonSkippingEvents = exonSkips.stream()
                .collect(Collectors.groupingBy(ExonSkip::getId, Collectors.counting()));
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(analysis))) {
            bw.write(output.toString());
            bw.write("# Gene statistics\n");
            bw.write("id\tname\tlength\texon skipping events\ttranscripts\tcds\texons\tintrons\n");
            for (Gene gene : gtfAnnotation.getGenes().values()) {
                long geneExonSkippingEevents = geneIdToTotalExonSkippingEvents.getOrDefault(gene.getGeneID(), 0L);
                bw.write(gene.getGeneID() + "\t" +
                        gene.getGeneName() + "\t" +
                        (gene.getInterval() != null ? gene.getInterval().getLength() : gene.getTranscripts().values().stream()
                                .flatMap(transcript -> transcript.getExons().stream())
                                .collect(Collectors.toCollection(TreeSet::new))
                                .stream()
                                .mapToInt(exon -> {
                                    int start = exon.getInterval().getStart();
                                    int end = exon.getInterval().getEnd();
                                    return end - start + 1;
                                })
                                .sum()) + "\t" +
                        geneExonSkippingEevents + "\t" +
                        gene.getTranscripts().size() + "\t" +
                        gene.getTranscripts().values().stream().mapToInt(transcript -> transcript.getCds().size()).sum() + "\t" +
                        gene.getTranscripts().values().stream().mapToInt(transcript -> transcript.getExons().size()).sum() + "\t" +
                        gene.getIntrons().size() + "\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ExonSkip exonSkip = (ExonSkip) o;
        return nprots == exonSkip.nprots && ntrans == exonSkip.ntrans && min_skipped_exon == exonSkip.min_skipped_exon && max_skipped_exon == exonSkip.max_skipped_exon && min_skipped_bases == exonSkip.min_skipped_bases && max_skipped_bases == exonSkip.max_skipped_bases && Objects.equals(id, exonSkip.id) && Objects.equals(symbol, exonSkip.symbol) && Objects.equals(chr, exonSkip.chr) && strand == exonSkip.strand && Objects.equals(SV, exonSkip.SV) && Objects.equals(WT, exonSkip.WT) && Objects.equals(SV_prots, exonSkip.SV_prots) && Objects.equals(WT_prots, exonSkip.WT_prots);
    }

    public String getChr() {
        return chr;
    }

    public String getId() {
        return id;
    }

    public int getMax_skipped_bases() {
        return max_skipped_bases;
    }

    public int getMax_skipped_exon() {
        return max_skipped_exon;
    }

    public int getMin_skipped_bases() {
        return min_skipped_bases;
    }

    public int getMin_skipped_exon() {
        return min_skipped_exon;
    }

    public int getNprots() {
        return nprots;
    }

    public int getNtrans() {
        return ntrans;
    }

    public StrandDirection getStrand() {
        return strand;
    }

    public Interval getSV() {
        return SV;
    }

    public Set<String> getSV_prots() {
        return SV_prots;
    }

    public String getSymbol() {
        return symbol;
    }

    public TreeSet<Interval> getWT() {
        return WT;
    }

    public Set<String> getWT_prots() {
        return WT_prots;
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, symbol, chr, strand, nprots, ntrans, SV, WT, SV_prots, WT_prots, min_skipped_exon, max_skipped_exon, min_skipped_bases, max_skipped_bases);
    }

    private record SkippedExonsBases(int minSkippedExons, int maxSkippedExons, int minSkippedBases,
                                     int maxSkippedBases) {
    }

    public static class Builder {
        private String id;
        private String symbol;
        private String chr;
        private StrandDirection strand;
        private int nprots;
        private int ntrans;
        private Interval SV;
        private TreeSet<Interval> WT;
        private Set<String> SV_prots = new HashSet<>();
        private Set<String> WT_prots = new HashSet<>();
        private int min_skipped_exon;
        private int max_skipped_exon;
        private int min_skipped_bases;
        private int max_skipped_bases;

        public Builder setId(String id) {
            this.id = id;
            return this;
        }

        public Builder setSymbol(String symbol) {
            this.symbol = symbol;
            return this;
        }

        public Builder setChr(String chr) {
            this.chr = chr;
            return this;
        }

        public Builder setStrand(StrandDirection strand) {
            this.strand = strand;
            return this;
        }

        public Builder setNprots(int nprots) {
            this.nprots = nprots;
            return this;
        }

        public Builder setNtrans(int ntrans) {
            this.ntrans = ntrans;
            return this;
        }

        public Builder setSV(Interval SV) {
            this.SV = SV;
            return this;
        }

        public Builder setWT(TreeSet<Interval> WT) {
            this.WT = WT;
            return this;
        }

        public Builder setSV_prots(Set<String> SV_prots) {
            this.SV_prots = SV_prots;
            return this;
        }

        public Builder setWT_prots(Set<String> WT_prots) {
            this.WT_prots = WT_prots;
            return this;
        }

        public Builder setMin_skipped_exon(int min_skipped_exon) {
            this.min_skipped_exon = min_skipped_exon;
            return this;
        }

        public Builder setMax_skipped_exon(int max_skipped_exon) {
            this.max_skipped_exon = max_skipped_exon;
            return this;
        }

        public Builder setMin_skipped_bases(int min_skipped_bases) {
            this.min_skipped_bases = min_skipped_bases;
            return this;
        }

        public Builder setMax_skipped_bases(int max_skipped_bases) {
            this.max_skipped_bases = max_skipped_bases;
            return this;
        }

        public ExonSkip build() {
            return new ExonSkip(this);
        }
    }
}
