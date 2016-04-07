/*
 *  MaskGenotypeTableBuilder
 * 
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.genotypecall.MaskGenotypeCallTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskGenotypeTableBuilder {

    private final byte DEFAULT_MASK_VALUE = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;

    private final GenotypeTable myBase;
    private final BitSet[] myBitSets;
    private byte myMaskValue = DEFAULT_MASK_VALUE;
    private final int myNumSites;
    private int myNextSite = 0;
    private int myNextTaxon = 0;
    private long myNextCount = 0;

    public MaskGenotypeTableBuilder(GenotypeTable base) {
        myBase = base;
        myNumSites = myBase.numberOfSites();
        int numTaxa = myBase.numberOfTaxa();
        myBitSets = new BitSet[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            myBitSets[t] = new OpenBitSet(myNumSites);
        }
    }

    public void setMaskValue(byte value) {
        myMaskValue = value;
    }

    public boolean get(int taxon, int site) {
        return myBitSets[taxon].fastGet(site);
    }

    public void set(int taxon, int site) {
        myBitSets[taxon].fastSet(site);
    }

    public void setNext(boolean value) {
        if (value) {
            myBitSets[myNextTaxon].fastSet(myNextSite);
        }
        myNextCount++;
        myNextTaxon = (int) (myNextCount / myNumSites);
        myNextSite = (int) (myNextCount % myNumSites);
    }

    public GenotypeTable build() {
        return new CoreGenotypeTable(
                new MaskGenotypeCallTable(myBase.genotypeMatrix(), myMaskValue, myBitSets),
                myBase.positions(),
                myBase.taxa(),
                myBase.depth(),
                myBase.alleleProbability(),
                myBase.referenceProbability(),
                myBase.dosage(),
                myBase.annotations()
        );
    }

}
