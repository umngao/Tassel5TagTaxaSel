/*
 *  ListSiteStats
 * 
 *  Created on Dec 5, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.lang.ref.WeakReference;
import java.util.AbstractList;
import java.util.HashMap;
import java.util.Map;
import net.maizegenetics.util.Tuple;

/**
 *
 * @author Terry Casstevens
 */
public class ListSiteStats extends AbstractList<Tuple<int[][], int[]>> {

    private final static Map<GenotypeCallTable, WeakReference<ListSiteStats>> myInstances = new HashMap<>();

    private final GenotypeCallTable myGenotype;
    private final Tuple<int[][], int[]>[] myCache;

    private ListSiteStats(GenotypeCallTable genotype) {
        myGenotype = genotype;
        myCache = new Tuple[myGenotype.numberOfSites()];
    }

    public static ListSiteStats getInstance(GenotypeCallTable genotype) {

        WeakReference<ListSiteStats> temp = myInstances.get(genotype);

        ListSiteStats result = null;
        if (temp != null) {
            result = temp.get();
        }

        if (result == null) {
            result = new ListSiteStats(genotype);
            myInstances.put(genotype, new WeakReference<>(result));
        }

        return result;

    }

    @Override
    public Tuple<int[][], int[]> get(int index) {
        if (myCache[index] == null) {
            myCache[index] = myGenotype.siteStats(index);
        }
        return myCache[index];
    }

    @Override
    public int size() {
        return myGenotype.numberOfSites();
    }

}
