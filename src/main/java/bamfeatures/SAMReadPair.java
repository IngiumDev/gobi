package bamfeatures;

import htsjdk.samtools.SAMRecord;
import readsimulator.IdenticalPair;

public class SAMReadPair extends IdenticalPair<SAMRecord> {


    public SAMReadPair(SAMRecord first, SAMRecord second) {
        super(first, second);
    }
}
