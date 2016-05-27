/*
 *  MaskSiteMatrix
 * 
 *  Created on May 6, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskSiteMatrix extends MaskMatrix {

    MaskSiteMatrix(BitSet[] bitSets, MaskMatrix merge) {
        super(bitSets, merge);
    }

    @Override
    public boolean get(int taxon, int site) {
        return myBitSets[site].fastGet(taxon);
    }

}
