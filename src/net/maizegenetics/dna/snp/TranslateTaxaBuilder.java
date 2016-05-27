/*
 *  TranslateTaxaBuilder
 * 
 *  Created on May 27, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateTaxaBuilder {

    private final BitSet myTaxaToKeep;
    private final TranslateTaxa myBase;
    private final int myNumBaseTaxa;

    public TranslateTaxaBuilder(int numBaseTaxa) {
        myNumBaseTaxa = numBaseTaxa;
        myTaxaToKeep = new OpenBitSet(numBaseTaxa);
        myBase = null;
    }

    public TranslateTaxaBuilder(TranslateTaxa base) {
        myNumBaseTaxa = base.numTaxa();
        myTaxaToKeep = new OpenBitSet(myNumBaseTaxa);
        myBase = base;
    }

    /**
     * Returns no translation instance.
     *
     * @param numTaxa number of taxa
     * @return no translation instance
     */
    public static TranslateTaxa noTranslation(int numTaxa) {
        return new TranslateTaxa(numTaxa);
    }

    public void keepTaxon(int taxon) {
        myTaxaToKeep.fastSet(taxon);
    }

    public TranslateTaxa build() {

        int numTaxaToKeep = (int) myTaxaToKeep.cardinality();

        if (numTaxaToKeep == 0) {
            throw new IllegalStateException("TranslateTaxaBuilder: build: no taxa to keep.");
        } else if (numTaxaToKeep == myNumBaseTaxa) {
            if (myBase != null) {
                return myBase;
            } else {
                return new TranslateTaxa(myNumBaseTaxa);
            }
        }

        int[] taxaRedirect = new int[numTaxaToKeep];
        int count = 0;
        if (myBase == null) {
            for (int i = 0; i < myNumBaseTaxa; i++) {
                if (myTaxaToKeep.fastGet(i)) {
                    taxaRedirect[count++] = i;
                }
            }
        } else {
            for (int i = 0; i < myNumBaseTaxa; i++) {
                if (myTaxaToKeep.fastGet(i)) {
                    taxaRedirect[count++] = myBase.translateTaxon(i);
                }
            }
        }

        return new TranslateTaxaRedirect(taxaRedirect);

    }

}
