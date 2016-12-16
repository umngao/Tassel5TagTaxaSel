/*
 *  TranslateSiteRedirect
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

import java.util.Arrays;

/**
 * Translation redirects site to corresponding site in base genotype table.
 *
 * @author Terry Casstevens
 */
public class TranslateSiteRedirect extends TranslateSite {

    protected final int[] mySiteRedirect;

    /**
     * Constructor
     *
     * @param siteRedirect redirect site indices to base site. Should be
     * ordered.
     */
    TranslateSiteRedirect(int[] siteRedirect) {
        super(siteRedirect.length);
        mySiteRedirect = siteRedirect;
    }

    /**
     * Translates site to base site.
     *
     * @param site site
     * @return translated base site
     */
    @Override
    public int translate(int site) {
        return mySiteRedirect[site];
    }

    /**
     * Translates base site to this site. Uses binary search algorithm since
     * indices are ordered.
     *
     * @param site site
     * @return translated site
     */
    @Override
    public int reverseTranslateSite(int site) {
        return Arrays.binarySearch(mySiteRedirect, site);
    }
    
    @Override
    public boolean hasTranslations() {
        return true;
    }

}
