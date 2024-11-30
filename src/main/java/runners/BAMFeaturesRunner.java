package runners;

import bamfeatures.ReadAnnotator;
import gtf.GTFTreeAnnotationTemp;
import gtf.structs.Gene;
import gtf.types.StrandDirection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.type.FileArgumentType;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import parsers.GTFParser;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

public class BAMFeaturesRunner {
    /* (BAM features) - Analyze mapped read annotation from BAM-files: (200pts):
Implement a program extracting information on all mapped read pairs.
Synopsis:
java -jar bamfeatures.jar -gtf <gtf_file> -bam <bamfile> -o <output_tsv>
[ -frstrand true/false ]
The option -frstrand specifies the strandness of the experiment. If the parameter is not given
the experiment is strand unspecific, if true the first read maps sense to the transcribed region, if
false the first read maps antisense to the transcribed region.
The output should contain one line per read pair in the format
readid<tab>annotation.
The annotation should consist of the tab separated features defined below. You find a table of
sample inputs along with the strandness information and path to the reference output in
/mnt/biosoft/praktikum/genprakt/BamFeatures/ref.table*/
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("BAMFeaturesRunner").build().defaultHelp(true)
                .description("Run BAMFeaturesRunner");
        parser.addArgument("-gtf").required(true).help("GTF file").metavar("<gtf_file>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-bam").required(true).help("BAM file").metavar("<bamfile>").type(new FileArgumentType().verifyIsFile());
        parser.addArgument("-o").required(true).help("Output file").metavar("<output_tsv>");
        parser.addArgument("-frstrand").help("true/false").metavar("<true/false>");
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

//        GTFTreeAnnotationTemp gtfTreeAnnotationTemp = new GTFTreeAnnotationTemp(GTFParser.parseGTF(res.getString("gtf")));
//        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(res.getString("bam")));
//        Iterator<SAMRecord> it = reader.iterator();
//        SAMRecord sr = it.next();
//        ArrayList<Gene> genes = gtfTreeAnnotationTemp.getTree("1", StrandDirection.REVERSE).getIntervalsIntersecting(24738, 24891, new ArrayList<>());
//        System.out.println();
        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(res.getString("bam")));
        StrandDirection strandSpecific = (res.getString("frstrand") == null) ? StrandDirection.UNSPECIFIED : (res.getString("frstrand").equals("true") ? StrandDirection.FORWARD : StrandDirection.REVERSE);
        ReadAnnotator readAnnotator = new ReadAnnotator.Builder()
                .setSamReader(reader)
                .setGtfFile(new File(res.getString("gtf")))
                .setOutputFile(new File(res.getString("o")))
                .setStrandSpecificity(strandSpecific)
                .build();
        readAnnotator.annotateReads();
        System.out.println();
    }
}
