/*
 *  TranslateSite
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

/**
 * No translation to site.
 * 
 * @author Terry Casstevens
 */
public class TranslateSite {

    private final int myNumSites;

    /**
     * Constructor
     *
     * @param numSites number of sites
     */
    TranslateSite(int numSites) {
        myNumSites = numSites;
    }

    /**
     * Translates site to base site. This class has no translation.
     *
     * @param site site
     * @return translated base site
     */
    public int translate(int site) {
        return site;
    }

    /**
     * Translates base site to this site. This class has no translation.
     *
     * @param site site
     * @return translated site
     */
    public int reverseTranslateSite(int site) {
        return site;
    }

    /**
     * Number of sites represented by this translation. Number of base sites
     * will be the same or larger.
     *
     * @return number of sites
     */
    public int numSites() {
        return myNumSites;
    }
    
    public boolean hasTranslations() {
        return false;
    }

}
