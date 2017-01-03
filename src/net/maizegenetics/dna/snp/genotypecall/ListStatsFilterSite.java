/*
 *  ListStatsFilterSite
 * 
 *  Created on Dec 27, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.Translate;
import net.maizegenetics.util.Tuple;

/**
 *
 * @author Terry Casstevens
 */
public class ListStatsFilterSite extends ListStats {

    private final ListStats myBase;
    private final Translate myTranslate;

    ListStatsFilterSite(FilterGenotypeCallTable genotype, ListStats base) {
        super(genotype, genotype.numberOfSites());
        myBase = base;
        myTranslate = genotype.myTranslate;
    }

    @Override
    public Tuple<int[][], int[]> get(int index) {
        Tuple<int[][], int[]> result = myBase.get(myTranslate.site(index));
        result.y[AlleleFreqCache.INDEX] = index;
        return result;
    }

}
