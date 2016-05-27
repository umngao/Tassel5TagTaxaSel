/*
 *  TranslateSiteBuilder
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
public class TranslateSiteBuilder {
    
    private final BitSet mySitesToKeep;
    private final TranslateSite myBase;
    private final int myNumBaseSites;
    
    public TranslateSiteBuilder(int numBaseSites) {
        myNumBaseSites = numBaseSites;
        mySitesToKeep = new OpenBitSet(numBaseSites);
        myBase = null;
    }
    
    public TranslateSiteBuilder(TranslateSite base) {
        myNumBaseSites = base.numSites();
        mySitesToKeep = new OpenBitSet(myNumBaseSites);
        myBase = base;
    }

    /**
     * Returns no translation instance.
     *
     * @param numSites number of sites
     * @return no translation instance
     */
    public static TranslateSite noTranslation(int numSites) {
        return new TranslateSite(numSites);
    }
    
    public static TranslateSite unorderedTranslation(int[] sitesNewOrder, TranslateSite base) {
        
        int numSites = sitesNewOrder.length;
        
        if (base == null) {
            int[] result = new int[numSites];
            System.arraycopy(sitesNewOrder, 0, result, 0, numSites);
            return new TranslateSiteRedirectUnordered(result);
        } else {
            if (numSites != base.numSites()) {
                throw new IllegalStateException("TranslateSiteBuilder: unorderedTranslation: number of newly ordered sites: " + numSites + " should equal base num sites: " + base.numSites());
            }
            int[] result = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = base.translateSite(sitesNewOrder[i]);
            }
            return new TranslateSiteRedirectUnordered(result);
        }
        
    }
    
    public void keepSite(int site) {
        mySitesToKeep.fastSet(site);
    }

    /**
     * Keeps a range of sites
     *
     * @param start lower index
     * @param end one-past the last site to keep
     */
    public void keepSites(int start, int end) {
        mySitesToKeep.set(start, end);
    }
    
    public TranslateSite build() {
        
        int numSitesToKeep = (int) mySitesToKeep.cardinality();
        
        if (numSitesToKeep == 0) {
            throw new IllegalStateException("TranslateSiteBuilder: build: no sites to keep.");
        } else if (numSitesToKeep == myNumBaseSites) {
            if (myBase != null) {
                return myBase;
            } else {
                return new TranslateSite(myNumBaseSites);
            }
        }
        
        int[] siteRedirect = new int[numSitesToKeep];
        int count = 0;
        if (myBase == null) {
            for (int i = 0; i < myNumBaseSites; i++) {
                if (mySitesToKeep.fastGet(i)) {
                    siteRedirect[count++] = i;
                }
            }
        } else {
            for (int i = 0; i < myNumBaseSites; i++) {
                if (mySitesToKeep.fastGet(i)) {
                    siteRedirect[count++] = myBase.translateSite(i);
                }
            }
        }
        
        return new TranslateSiteRedirect(siteRedirect);
        
    }
    
}
