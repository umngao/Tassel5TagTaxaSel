/*
 *  SiteScoreBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Writer;

import java.util.LinkedHashMap;
import java.util.Map;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;

/**
 *
 * @author Terry Casstevens
 */
public class SiteScoreBuilder {

    private final Map<GenotypeTable.SITE_SCORE_TYPE, Byte2DBuilder> myBuilders = new LinkedHashMap<>();
    private final boolean myIsHDF5;

    private SiteScoreBuilder(GenotypeTable.SITE_SCORE_TYPE[] siteScoreTypes, int numTaxa, int numSites) {
        for (int i = 0; i < siteScoreTypes.length; i++) {
            myBuilders.put(siteScoreTypes[i], Byte2DBuilder.getInstance(numTaxa, numSites, siteScoreTypes[i]));
        }
        myIsHDF5 = false;
    }

    private SiteScoreBuilder(IHDF5Writer writer, GenotypeTable.SITE_SCORE_TYPE[] siteScoreTypes, int numTaxa, int numSites) {
        for (int i = 0; i < siteScoreTypes.length; i++) {
            myBuilders.put(siteScoreTypes[i], Byte2DBuilder.getInstance(writer, numSites, siteScoreTypes[i]));
        }
        myIsHDF5 = false;
    }

    public SiteScoreBuilder getInstance(GenotypeTable.SITE_SCORE_TYPE[] siteScoreTypes, int numTaxa, int numSites) {
        return new SiteScoreBuilder(siteScoreTypes, numTaxa, numSites);
    }

    public SiteScoreBuilder getInstance(IHDF5Writer writer, GenotypeTable.SITE_SCORE_TYPE[] siteScoreTypes, int numTaxa, int numSites) {
        return null;
    }

    public SiteScoreBuilder set(GenotypeTable.SITE_SCORE_TYPE siteScoreType, int taxon, int site, float value) {
        if (myIsHDF5) {
            throw new IllegalStateException("SiteScoreBuilder: set: use addTaxon for HDF5 files.");
        }
        myBuilders.get(siteScoreType).set(taxon, site, SiteScoreUtil.floatToBytePercentage(value));
        return this;
    }

    public SiteScore build() {
        Byte2D[] input = new Byte2D[myBuilders.size()];
        int count = 0;
        for (Byte2DBuilder builder : myBuilders.values()) {
            input[count++] = builder.build();
        }
        myBuilders.clear();
        return new SiteScore(input);
    }

}
