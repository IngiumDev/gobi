package me.ingiumdev.gobi.runners;

import me.ingiumdev.gobi.go.GeneSetEnrichmentAnalysis;
import me.ingiumdev.gobi.go.RootType;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.File;

import static me.ingiumdev.gobi.go.RootType.parseRootType;

public class GOEnrichmentRunner {
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("GOEnrichmentRunner").build().defaultHelp(true)
                .description("Perform Gene Ontology Set Enrichment Analysis");

        parser.addArgument("-obo")
                .required(true)
                .help("Path to the OBO file")
                .metavar("<obo_file>")
                .type(new FileArgumentType().verifyIsFile());

        parser.addArgument("-root")
                .required(true)
                .help("Specify GO namespace: biological_process, cellular_component, or molecular_function")
                .metavar("<GO_namespace>")
                .choices("biological_process", "cellular_component", "molecular_function");

        parser.addArgument("-mapping")
                .required(true)
                .help("Path to the gene-to-GO mapping file")
                .metavar("<gene2go_mapping_file>")
                .type(new FileArgumentType().verifyIsFile());

        parser.addArgument("-mappingtype")
                .required(true)
                .choices("ensembl", "go")
                .help("Specify the format of the mapping file: ensembl or go")
                .metavar("[ensembl|go]");

        parser.addArgument("-enrich")
                .required(true)
                .help("Path to the enrichment analysis input file")
                .metavar("<diffexp_file>")
                .type(new FileArgumentType().verifyIsFile());

        parser.addArgument("-o")
                .required(true)
                .help("Path to the output TSV file for enrichment results")
                .metavar("<output_tsv>");

        parser.addArgument("-minsize")
                .required(true)
                .help("Minimum size of GO entries to consider")
                .metavar("<int>")
                .type(Integer.class);

        parser.addArgument("-maxsize")
                .required(true)
                .help("Maximum size of GO entries to consider")
                .metavar("<int>")
                .type(Integer.class);

        parser.addArgument("-overlapout")
                .help("Optional output file for DAG overlap information")
                .metavar("<overlap_out_tsv>");

        if (args.length == 0) {
            parser.printHelp();
            System.exit(1);
        }

        try {
            Namespace res = parser.parseArgs(args);
            if ((Integer) res.get("minsize") > (Integer) res.get("maxsize")) {
                System.err.println("Error: minsize must be less than or equal to maxsize.");
                System.exit(1);
            }
            start(res);
        } catch (ArgumentParserException e) {
            parser.printHelp();
            System.exit(1);
        }
    }

    private static void start(Namespace res) {
        // Extract parsed arguments
        String oboFile = res.getString("obo");
        String mappingFile = res.getString("mapping");
        RootType root = parseRootType(res.get("root"));
        String mappingType = res.get("mappingtype");
        String enrichFile = res.getString("enrich");
        File outputTsv = new File(res.getString("o"));
        int minSize = res.get("minsize");
        int maxSize = res.get("maxsize");
        String overlapOut = res.get("overlapout");
        GeneSetEnrichmentAnalysis analysis = new GeneSetEnrichmentAnalysis.Builder().setRootType(root)
                .setMinSize(minSize)
                .setMaxSize(maxSize)
                .build();
        analysis.initMapping(mappingFile);
        analysis.initDAG(oboFile);
        analysis.initDifferentialExpression(enrichFile);
    }
}
