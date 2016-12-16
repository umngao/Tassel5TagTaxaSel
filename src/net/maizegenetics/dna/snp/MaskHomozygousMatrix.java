/*
 *  MaskHomozygousMatrix
 * 
 *  Created on Dec 12, 2016
 */
package net.maizegenetics.dna.snp;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.ForkJoinPool;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskHomozygousMatrix implements MaskMatrix {

    private final GenotypeCallTable myGenotype;
    private final int myNumTaxa;
    private final int myNumSites;
    private final Cache<Integer, BitSet> myCache;
    private final CopyOnWriteArraySet<Integer> myCurrentlyProcessingBlocks = new CopyOnWriteArraySet<>();
    private final ForkJoinPool myThreadPool = new ForkJoinPool();

    MaskHomozygousMatrix(GenotypeCallTable genotype) {
        myGenotype = genotype;
        myNumTaxa = genotype.numberOfTaxa();
        myNumSites = genotype.numberOfSites();
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
            byte[] temp = myGenotype.genotypeForAllTaxa(site);
            for (int t = 0; t < myNumTaxa; t++) {
                if ((temp[t] & 0xF) != (temp[t] >>> 4)) {
                    result.fastSet(t);
                }
            }
            myCache.put(site, result);
            return result;
        }
    }

    @Override
    public boolean get(int taxon, int site) {
        
        BitSet temp = myCache.getIfPresent(site);
        if (temp == null) {
            if (myCurrentlyProcessingBlocks.add(site)) {
                myThreadPool.submit(new ProcessSite(site));
            }
        } else {
            return temp.fastGet(taxon);
        }
        
        byte value = myGenotype.genotype(taxon, site);
        return (value & 0xF) != (value >>> 4);
        
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
            if ((temp[s] & 0xF) != (temp[s] >>> 4)) {
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

    private class ProcessSite implements Runnable {

        private final int mySite;

        public ProcessSite(int site) {
            mySite = site;
        }

        @Override
        public void run() {
            BitSet result = new OpenBitSet(myNumTaxa);
            byte[] temp = myGenotype.genotypeForAllTaxa(mySite);
            for (int t = 0; t < myNumTaxa; t++) {
                if ((temp[t] & 0xF) != (temp[t] >>> 4)) {
                    result.fastSet(t);
                }
            }
            myCache.put(mySite, result);
            myCurrentlyProcessingBlocks.remove(mySite);
        }

    }

}
