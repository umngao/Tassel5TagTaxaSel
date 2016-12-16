/*
 *  TranslateSiteBuilder
 * 
 *  Created on May 27, 2016
 */
package net.maizegenetics.dna.snp;

import java.util.Arrays;
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

    private TranslateSiteBuilder(int numBaseSites) {
        myNumBaseSites = numBaseSites;
        mySitesToKeep = new OpenBitSet(numBaseSites);
        myBase = null;
    }

    private TranslateSiteBuilder(TranslateSite base) {
        myNumBaseSites = base.numSites();
        mySitesToKeep = new OpenBitSet(myNumBaseSites);
        myBase = base;
    }

    public static TranslateSiteBuilder getInstance(int numBaseSites) {
        return new TranslateSiteBuilder(numBaseSites);
    }

    public static TranslateSiteBuilder getInstance(TranslateSite base) {
        return new TranslateSiteBuilder(base);
    }

    public static TranslateSiteBuilder getInstance(int numBaseSites, TranslateSite base) {
        if (base == null) {
            return new TranslateSiteBuilder(numBaseSites);
        } else {
            if (numBaseSites != base.numSites()) {
                throw new IllegalArgumentException("TranslateTaxaBuilder: getInstance: numBaseSites: " + numBaseSites + " should equal base: " + base.numSites());
            }
            return new TranslateSiteBuilder(base);
        }
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

    public static TranslateSite orderedTranslation(int[] siteRedirect, TranslateSite base) {

        int numSites = siteRedirect.length;

        if (base == null) {
            int[] result = new int[numSites];
            System.arraycopy(siteRedirect, 0, result, 0, numSites);
            Arrays.sort(result);
            return new TranslateSiteRedirectUnordered(result);
        } else {
            int[] result = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = base.translate(siteRedirect[i]);
            }
            Arrays.sort(result);
            return new TranslateSiteRedirectUnordered(result);
        }

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
                result[i] = base.translate(sitesNewOrder[i]);
            }
            return new TranslateSiteRedirectUnordered(result);
        }

    }

    /**
     * Keeps a range of sites from start (inclusive) to end (inclusive)
     * 
     * @param start start site
     * @param end end site
     * @param base base translation
     * 
     * @return new translation
     */
    public static TranslateSite range(int start, int end, TranslateSite base) {

        if (base == null) {
            return new TranslateSiteRange(start, end);
        } else {
            return new TranslateSiteBuilder(base).keepSites(start, end).build();
        }

    }

    /**
     * Keep specified site.
     *
     * @param site site
     *
     * @return this builder
     */
    public TranslateSiteBuilder keepSite(int site) {
        mySitesToKeep.fastSet(site);
        return this;
    }

    /**
     * Keeps a range of sites from start (inclusive) to end (inclusive)
     *
     * @param start lower index
     * @param end last site to keep
     *
     * @return this builder
     */
    public TranslateSiteBuilder keepSites(int start, int end) {
        mySitesToKeep.set(start, end + 1);
        return this;
    }

    public TranslateSiteBuilder keepSites(int[] sites) {
        for (int current : sites) {
            keepSite(current);
        }
        return this;
    }

    public int numSites() {
        return (int) mySitesToKeep.cardinality();
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
                    siteRedirect[count++] = myBase.translate(i);
                }
            }
        }

        // Sites in siteRedirect already sorted since they were
        // put into the array in order.
        return new TranslateSiteRedirect(siteRedirect);

    }

}
