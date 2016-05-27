/*
 *  TranslateSiteRedirectUnordered
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateSiteRedirectUnordered extends TranslateSiteRedirect {

    /**
     * Constructor
     *
     * @param siteRedirect redirected site indices to base site. These are
     * unordered.
     */
    TranslateSiteRedirectUnordered(int[] siteRedirect) {
        super(siteRedirect);
    }

    /**
     * Translates base site to this site. Searches all indices since not
     * ordered.
     *
     * @param site site
     * @return translated site
     */
    @Override
    public int reverseTranslateSite(int site) {
        for (int i = 0; i < numSites(); i++) {
            if (mySiteRedirect[i] == site) {
                return i;
            }
        }
        return -1;
    }

}
