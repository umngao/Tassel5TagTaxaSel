/*
 *  AlleleDepthBuilder
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
public class AlleleDepthBuilder {

    private final Map<SiteScore.SITE_SCORE_TYPE, Byte2DBuilder> myBuilders = new LinkedHashMap<>();
    private final int myNumSites;

    private AlleleDepthBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < AlleleDepth.ALLELE_DEPTH_TYPES.length; i++) {
            myBuilders.put(AlleleDepth.ALLELE_DEPTH_TYPES[i], Byte2DBuilder.getInstance(numTaxa, numSites, AlleleDepth.ALLELE_DEPTH_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    private AlleleDepthBuilder(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        for (int i = 0; i < AlleleDepth.ALLELE_DEPTH_TYPES.length; i++) {
            myBuilders.put(AlleleDepth.ALLELE_DEPTH_TYPES[i], Byte2DBuilder.getInstance(writer, numSites, AlleleDepth.ALLELE_DEPTH_TYPES[i], taxaList));
        }
        myNumSites = numSites;
    }

    public static AlleleDepthBuilder getInstance(IHDF5Writer writer, int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleDepthBuilder(writer, numTaxa, numSites, taxaList);
    }

    public static AlleleDepthBuilder getAlleleDepthInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new AlleleDepthBuilder(numTaxa, numSites, taxaList);
    }

    public static AlleleDepth getFilteredInstance(AlleleDepth base, FilterGenotypeTable filterGenotypeTable) {
        Collection<Byte2D> storage = base.byteStorage();
        FilterByte2D[] resultStorage = new FilterByte2D[storage.size()];
        int count = 0;
        for (Byte2D current : storage) {
            resultStorage[count] = Byte2DBuilder.getFilteredInstance(current, filterGenotypeTable);
        }
        return new AlleleDepth(resultStorage);
    }

    public AlleleDepthBuilder addTaxon(int taxon, int[] values, SiteScore.SITE_SCORE_TYPE type) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("AlleleDepthBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        byte[] result = AlleleDepthUtil.depthIntToByte(values);
        myBuilders.get(type).addTaxon(taxon, result);
        return this;
    }

    public AlleleDepth build() {
        Byte2D[] input = new Byte2D[myBuilders.size()];
        int count = 0;
        for (Byte2DBuilder builder : myBuilders.values()) {
            input[count++] = builder.build();
        }
        myBuilders.clear();
        return new AlleleDepth(input);
    }

}
