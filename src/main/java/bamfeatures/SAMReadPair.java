package bamfeatures;

import gtf.structs.Interval;
import htsjdk.samtools.SAMRecord;
import readsimulator.IdenticalPair;

import java.util.TreeSet;

public class SAMReadPair extends IdenticalPair<SAMRecord> {


    public SAMReadPair(SAMRecord first, SAMRecord second) {
        super(first, second);
    }
}
