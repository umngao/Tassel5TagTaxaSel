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
import net.maizegenetics.dna.snp.GenotypeTable;
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

    public byte majorAllele(int site) {
        int[][] alleles = get(site).x;
        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public double majorAlleleFrequency(int site) {

        int[][] alleles = get(site).x;
        int numAlleles = alleles[0].length;
        if (numAlleles >= 1) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing += alleles[1][i];
            }
            return (double) alleles[1][0] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    public byte minorAllele(int site) {
        int[][] alleles = get(site).x;
        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public double minorAlleleFrequency(int site) {

        int[][] alleles = get(site).x;
        int numAlleles = alleles[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing += alleles[1][i];
            }
            return (double) alleles[1][1] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

}
