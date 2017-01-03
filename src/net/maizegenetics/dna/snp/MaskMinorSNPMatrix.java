/*
 *  MaskMinorSNPMatrix
 * 
 *  Created on Dec 12, 2016
 */
package net.maizegenetics.dna.snp;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.ListStats;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskMinorSNPMatrix implements MaskMatrix {

    private final GenotypeCallTable myGenotype;
    private final int myNumTaxa;
    private final int myNumSites;
    private final Cache<Integer, BitSet> myCache;
    private final ListStats myStats;

    MaskMinorSNPMatrix(GenotypeCallTable genotype) {
        myGenotype = genotype;
        myNumTaxa = genotype.numberOfTaxa();
        myNumSites = genotype.numberOfSites();
        myStats = ListStats.getSiteInstance(genotype);
        myCache = CacheBuilder.newBuilder()
                .initialCapacity(1000)
                .maximumSize(1000)
                .build();
    }

    private BitSet getFromCache(int site) {
        BitSet result = myCache.getIfPresent(site);
        if (result != null) {
            return result;
        } else {
            result = new OpenBitSet(myNumTaxa);
            byte major = myStats.majorAllele(site);
            byte minor = myStats.minorAllele(site);
            byte[] temp = myGenotype.genotypeForAllTaxa(site);
            for (int t = 0; t < myNumTaxa; t++) {
                if (((temp[t] & 0xF) != major) && ((temp[t] & 0xF) != minor)) {
                    result.fastSet(t);
                } else if (((temp[t] >>> 4) != major) && ((temp[t] >>> 4) != minor)) {
                    result.fastSet(t);
                }
            }
            myCache.put(site, result);
            return result;
        }
    }

    @Override
    public boolean get(int taxon, int site) {
        return getFromCache(site).fastGet(taxon);
    }

    @Override
    public boolean isTaxonMaskedHint(int taxon) {
        return true;
    }

    @Override
    public BitSet maskForTaxon(int taxon) {
        BitSet result = new OpenBitSet(myNumSites);
        byte[] temp = myGenotype.genotypeForAllSites(taxon);
        for (int s = 0; s < myNumSites; s++) {
            byte major = myStats.majorAllele(s);
            byte minor = myStats.minorAllele(s);
            if (((temp[s] & 0xF) != major) && ((temp[s] & 0xF) != minor)) {
                result.fastSet(s);
            } else if (((temp[s] >>> 4) != major) && ((temp[s] >>> 4) != minor)) {
                result.fastSet(s);
            }
        }
        return result;
    }

    @Override
    public boolean isSiteMaskedHint(int site) {
        return getFromCache(site).cardinality() != 0;
    }

    @Override
    public BitSet maskForSite(int site) {
        return UnmodifiableBitSet.getInstance(getFromCache(site));
    }

    @Override
    public int numTaxa() {
        return myNumTaxa;
    }

    @Override
    public int numSites() {
        return myNumSites;
    }

    @Override
    public boolean isSiteOptimized() {
        return true;
    }

}
