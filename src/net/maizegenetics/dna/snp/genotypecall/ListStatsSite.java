/*
 *  ListStatsSite
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.util.Tuple;

/**
 *
 * @author Terry Casstevens
 */
public class ListStatsSite extends ListStats {

    private final Tuple<int[][], int[]>[] myCache;

    ListStatsSite(GenotypeCallTable genotype) {
        super(genotype, genotype.numberOfSites());
        myCache = new Tuple[genotype.numberOfSites()];
    }

    @Override
    public Tuple<int[][], int[]> get(int index) {
        if (myCache[index] == null) {
            myCache[index] = myGenotype.siteStats(index);
        }
        int num = myCache[index].x[0].length;
        int[][] alleleCounts = new int[2][num];
        System.arraycopy(myCache[index].x[0], 0, alleleCounts[0], 0, num);
        System.arraycopy(myCache[index].x[1], 0, alleleCounts[1], 0, num);
        num = myCache[index].y.length;
        int[] stats = new int[num];
        System.arraycopy(myCache[index].y, 0, stats, 0, num);
        return new Tuple<>(alleleCounts, stats);
    }

}
