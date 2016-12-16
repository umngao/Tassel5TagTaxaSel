/*
 *  TranslateTaxaRedirect
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateTaxaRedirect extends TranslateTaxa {

    private final int[] myTaxaRedirect;

    /**
     * Constructor
     *
     * @param taxaRedirect taxa redirected indices to base taxa. Should be
     * ordered.
     */
    TranslateTaxaRedirect(int[] taxaRedirect) {
        super(taxaRedirect.length);
        myTaxaRedirect = taxaRedirect;
    }

    /**
     * Returns translated taxon.
     *
     * @param taxon taxon
     * @return translated taxon
     */
    @Override
    public int translate(int taxon) {
        return myTaxaRedirect[taxon];
    }
    
    @Override
    public boolean hasTranslations() {
        return true;
    }

}
