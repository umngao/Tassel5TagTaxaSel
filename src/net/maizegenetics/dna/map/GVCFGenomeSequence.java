package net.maizegenetics.dna.map;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by zrm22 on 3/27/17.
 *
 * Interface is used to store GenomeSequences defined by a GATK generated GVCF file
 */
public interface GVCFGenomeSequence extends GenomeSequence {
    public HashMap<Chromosome,ArrayList<ArrayList<Integer>>> getConsecutiveRegions();
    public void writeFASTA(String fileName);
}
