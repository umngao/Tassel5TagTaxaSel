/*
 *  TranslateTaxa
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateTaxa {

    private final int myNumTaxa;

    /**
     * Constructor
     *
     * @param numTaxa number of taxa
     */
    TranslateTaxa(int numTaxa) {
        myNumTaxa = numTaxa;
    }

    /**
     * Translates taxon to base taxon. This class has no translation.
     *
     * @param taxon taxon
     * @return translated taxon
     */
    public int translate(int taxon) {
        return taxon;
    }

    /**
     * Number of taxa represented by this translation. Number of base taxa will
     * be the same or larger.
     *
     * @return number of taxa
     */
    public int numTaxa() {
        return myNumTaxa;
    }
    
    public boolean hasTranslations() {
        return false;
    }

    public TaxaList taxa(TaxaList orig) {
        TaxaListBuilder builder = new TaxaListBuilder();
        for (int t = 0; t < myNumTaxa; t++) {
            int current = translate(t);
            if (current == -1) {
                builder.add(new Taxon("Taxon" + t));
            } else {
                builder.add(orig.get(current));
            }
        }
        return builder.build();
    }

}
