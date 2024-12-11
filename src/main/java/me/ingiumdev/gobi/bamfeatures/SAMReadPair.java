package me.ingiumdev.gobi.bamfeatures;

import htsjdk.samtools.SAMRecord;
import me.ingiumdev.gobi.readsimulator.IdenticalPair;

public class SAMReadPair extends IdenticalPair<SAMRecord> {


    public SAMReadPair(SAMRecord first, SAMRecord second) {
        super(first, second);
    }
}
