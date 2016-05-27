/*
 *  MaskMatrix
 * 
 *  Created on May 6, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskMatrix {

    protected final BitSet[] myBitSets;

    /**
     * Constructs a MaskMatrix for use with components of a GenotypeTable.
     *
     * @param bitSets set bits to indicate which are masked
     * @param merge merge mask to prevent nested masks. can be null.
     */
    MaskMatrix(BitSet[] bitSets, MaskMatrix merge) {
        myBitSets = bitSets;
        if (merge != null) {
            int num = myBitSets.length;
            for (int i = 0; i < num; i++) {
                myBitSets[i].or(merge.myBitSets[i]);
            }
        }
    }

    public boolean get(int taxon, int site) {
        return myBitSets[taxon].fastGet(site);
    }

}
