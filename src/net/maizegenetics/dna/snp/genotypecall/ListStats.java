/*
 *  ListStats
 * 
 *  Created on Dec 27, 2016
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
public abstract class ListStats extends AbstractList<Tuple<int[][], int[]>> {

    private final static Map<GenotypeCallTable, WeakReference<ListStats>> TAXA_INSTANCES = new HashMap<>();
    private final static Map<GenotypeCallTable, WeakReference<ListStats>> SITE_INSTANCES = new HashMap<>();

    protected final GenotypeCallTable myGenotype;
    private final int myNumIndices;

    ListStats(GenotypeCallTable genotype, int numIndices) {
        myGenotype = genotype;
        myNumIndices = numIndices;
    }

    public static ListStats getTaxaInstance(GenotypeCallTable genotype) {

        WeakReference<ListStats> temp = TAXA_INSTANCES.get(genotype);

        ListStats result = null;
        if (temp != null) {
            result = temp.get();
        }

        if (result == null) {
            if (genotype instanceof FilterGenotypeCallTable) {
                FilterGenotypeCallTable filter = (FilterGenotypeCallTable) genotype;
                WeakReference<ListStats> temp1 = TAXA_INSTANCES.get(filter.myBaseGenotype);
                ListStats base = null;
                if (temp1 != null) {
                    base = temp1.get();
                }
                if (base == null || filter.myTranslate.hasSiteTranslations()) {
                    result = new ListStatsTaxa(genotype);
                } else {
                    result = new ListStatsFilterSite(filter, base);
                }
            } else {
                result = new ListStatsTaxa(genotype);
            }
            TAXA_INSTANCES.put(genotype, new WeakReference<>(result));
        }

        return result;

    }

    public static ListStats getSiteInstance(GenotypeCallTable genotype) {

        WeakReference<ListStats> temp = SITE_INSTANCES.get(genotype);

        ListStats result = null;
        if (temp != null) {
            result = temp.get();
        }

        if (result == null) {
            System.out.println("result is null");
            if (genotype instanceof FilterGenotypeCallTable) {
                FilterGenotypeCallTable filter = (FilterGenotypeCallTable) genotype;
                WeakReference<ListStats> temp1 = SITE_INSTANCES.get(filter.myBaseGenotype);
                ListStats base = null;
                if (temp1 != null) {
                    base = temp1.get();
                }
                if (base == null) {
                System.out.println("base is null");
                }
                if (filter.myTranslate.hasTaxaTranslations()) {
                    System.out.println("has taxa translations");
                }
                if (base == null || filter.myTranslate.hasTaxaTranslations()) {
                    result = new ListStatsSite(genotype);
                    System.out.println("result is ListStatsSite with Filter");
                } else {
                    result = new ListStatsFilterSite(filter, base);
                    System.out.println("result is ListStatsFilterSite");
                }
            } else {
                result = new ListStatsSite(genotype);
                System.out.println("result is ListStatsSite");
            }
            SITE_INSTANCES.put(genotype, new WeakReference<>(result));
        }

        return result;

    }

    @Override
    public int size() {
        return myNumIndices;
    }

    public byte majorAllele(int index) {
        int[][] alleles = get(index).x;
        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public double majorAlleleFrequency(int index) {

        int[][] alleles = get(index).x;
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

    public byte minorAllele(int index) {
        int[][] alleles = get(index).x;
        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public double minorAlleleFrequency(int index) {

        int[][] alleles = get(index).x;
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
