/*
 *  MaskGenotypeCallTable
 * 
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.util.BitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskGenotypeCallTable extends AbstractGenotypeCallTable {
    
    private final GenotypeCallTable myBase;
    private final byte myMaskValue;
    private final BitSet[] myBitSets;

    public MaskGenotypeCallTable(GenotypeCallTable base, byte maskValue, BitSet[] bitSets) {
        super(base.numberOfTaxa(), base.numberOfSites(), base.isPhased(), base.alleleDefinitions());
        myBase = base;
        myMaskValue = maskValue;
        myBitSets = bitSets;
    }

    @Override
    public byte genotype(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBase.transposeData(siteInnerLoop);
    }

}
