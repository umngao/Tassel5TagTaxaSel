/*
 *  ListStatsFilterTaxa
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
public class ListStatsFilterTaxa extends ListStats {

    private final ListStats myBase;
    private final Translate myTranslate;

    ListStatsFilterTaxa(FilterGenotypeCallTable genotype, ListStats base) {
        super(genotype, genotype.numberOfTaxa());
        myBase = base;
        myTranslate = genotype.myTranslate;
    }

    @Override
    public Tuple<int[][], int[]> get(int index) {
        Tuple<int[][], int[]> result = myBase.get(myTranslate.taxon(index));
        result.y[AlleleFreqCache.INDEX] = index;
        return result;
    }

}
