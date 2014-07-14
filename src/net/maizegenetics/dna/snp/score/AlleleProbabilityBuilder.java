/*
 *  AlleleProbabilityBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Writer;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.dna.snp.byte2d.FilterByte2D;
import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleProbabilityBuilder {

    private final Map<SiteScore.SITE_SCORE_TYPE, Byte2DBuilder> myBuilders = new LinkedHashMap<>();
    private final int myNumSites;

    private static final SiteScore.SITE_SCORE_TYPE[] ALLELE_PROBABILITY_TYPES = new SiteScore.SITE_SCORE_TYPE[]{
        SiteScore.SITE_SCORE_TYPE.ProbA, SiteScore.SITE_SCORE_TYPE.ProbC,
        SiteScore.SITE_SCORE_TYPE.ProbG, SiteScore.SITE_SCORE_TYPE.ProbT,
        SiteScore.SITE_SCORE_TYPE.ProbGap, SiteScore.SITE_SCORE_TYPE.ProbInsertion};

    private AlleleProbabilityBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < ALLELE_PROBABILITY_TYPES.length; i++) {
            myBuilders.put(ALLELE_PROBABILITY_TYPES[i], Byte2DBuilder.getInstance(numTaxa, numSites, ALLELE_PROBABILITY_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    private AlleleProbabilityBuilder(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < ALLELE_PROBABILITY_TYPES.length; i++) {
            myBuilders.put(ALLELE_PROBABILITY_TYPES[i], Byte2DBuilder.getInstance(writer, numSites, ALLELE_PROBABILITY_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    public static AlleleProbabilityBuilder getInstance(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleProbabilityBuilder(writer, numTaxa, numSites, taxaList);
    }

    public static AlleleProbabilityBuilder getAlleleProbabilityInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleProbabilityBuilder(numTaxa, numSites, taxaList);
    }

    public static AlleleProbability getFilteredInstance(AlleleProbability base, FilterGenotypeTable filterGenotypeTable) {
        Collection<Byte2D> storage = base.byteStorage();
        FilterByte2D[] resultStorage = new FilterByte2D[storage.size()];
        int count = 0;
        for (Byte2D current : storage) {
            resultStorage[count] = Byte2DBuilder.getFilteredInstance(current, filterGenotypeTable);
        }
        return new AlleleProbability(resultStorage);
    }

    public AlleleProbabilityBuilder addTaxon(int taxon, byte[] values, SiteScore.SITE_SCORE_TYPE type) {
        myBuilders.get(type).addTaxon(taxon, values);
        return this;
    }

    public AlleleProbabilityBuilder addTaxon(int taxon, float[] values, SiteScore.SITE_SCORE_TYPE type) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("AlleleProbabilityBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        byte[] result = SiteScoreUtil.floatToBytePercentage(values);
        myBuilders.get(type).addTaxon(taxon, result);
        return this;
    }

    public AlleleProbability build() {
        Byte2D[] input = new Byte2D[myBuilders.size()];
        int count = 0;
        for (Byte2DBuilder builder : myBuilders.values()) {
            input[count++] = builder.build();
        }
        myBuilders.clear();
        return new AlleleProbability(input);
    }

}
