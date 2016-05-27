/*
 *  TranslateTaxa
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

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
    public int translateTaxon(int taxon) {
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

}
