/*
 *  ImputeProbabilityBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Writer;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.dna.snp.byte2d.FilterByte2D;
import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class ImputeProbabilityBuilder {

    private Byte2DBuilder myBuilder;

    private ImputeProbabilityBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(numTaxa, numSites, GenotypeTable.SITE_SCORE_TYPE.ImputedProbablity, taxaList);
    }

    private ImputeProbabilityBuilder(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(writer, numSites, GenotypeTable.SITE_SCORE_TYPE.ImputedProbablity, taxaList);
    }

    public static ImputeProbabilityBuilder getInstance(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        return new ImputeProbabilityBuilder(writer, numTaxa, numSites, taxaList);
    }

    public static ImputeProbabilityBuilder getInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new ImputeProbabilityBuilder(numTaxa, numSites, taxaList);
    }

    public static ImputeProbability getFilteredInstance(ImputeProbability base, FilterGenotypeTable filterGenotypeTable) {
        FilterByte2D resultStorage = Byte2DBuilder.getFilteredInstance(base.byteStorage(GenotypeTable.SITE_SCORE_TYPE.ImputedProbablity), filterGenotypeTable);
        return new ImputeProbability(resultStorage);
    }

    public ImputeProbabilityBuilder addTaxon(int taxon, byte[] values, GenotypeTable.SITE_SCORE_TYPE type) {
        myBuilder.addTaxon(taxon, values);
        return this;
    }

    public ImputeProbability build() {
        Byte2D input = myBuilder.build();
        myBuilder = null;
        return new ImputeProbability(input);
    }

}
