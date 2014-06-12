/*
 *  SiteScore
 */
package net.maizegenetics.dna.snp.score;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class SiteScore {

    private final Map<GenotypeTable.SITE_SCORE_TYPE, Byte2D> myValues = new HashMap<>();
    private final GenotypeTable.SITE_SCORE_TYPE myOnlyScoreType;
    private final int myNumTaxa;
    private final int myNumSites;

    SiteScore(Byte2D[] values) {
        if (values.length == 0) {
            throw new IllegalArgumentException("SiteScore: init: no values provided.");
        }
        int numTaxa = values[0].numTaxa();
        int numSites = values[0].numSites();
        for (int i = 0; i < values.length; i++) {
            if ((numTaxa != values[i].numTaxa()) || (numSites != values[i].numSites())) {
                throw new IllegalArgumentException("SiteScore: init: number of taxa or sites don't match for all values.");
            }
            myValues.put(values[i].siteScoreType(), values[i]);
        }
        myNumTaxa = numTaxa;
        myNumSites = numSites;
        if (values.length == 1) {
            myOnlyScoreType = values[0].siteScoreType();
        } else {
            myOnlyScoreType = null;
        }
    }

    public float siteScore(int taxon, int site, GenotypeTable.SITE_SCORE_TYPE scoreType) {
        return SiteScoreUtil.byteToFloatPercentage(myValues.get(scoreType).valueForAllele(taxon, site));
    }

    /**
     * Returns the site score of the given taxon and site. Used if only one site
     * score type represented by this instance.
     *
     * @param taxon taxon index
     * @param site site
     *
     * @return site score.
     */
    public float siteScore(int taxon, int site) {
        return siteScore(taxon, site, myOnlyScoreType);
    }

    /**
     * Returns the site scores (First dimension is taxa and second dimension is
     * sites. Default is allele 0. Can be used with scores that only have one
     * value per taxon and site.
     *
     * @return site scores.
     */
    public float[][] siteScores() {
        float[][] result = new float[numTaxa()][numSites()];
        for (int t = 0; t < numTaxa(); t++) {
            for (int s = 0; s < numSites(); s++) {
                result[t][s] = siteScore(t, s, myOnlyScoreType);
            }
        }
        return result;
    }

    /**
     * Return the site scores types.
     *
     * @return site score types.
     */
    public Set<GenotypeTable.SITE_SCORE_TYPE> siteScoreTypes() {
        return myValues.keySet();
    }

    public int numTaxa() {
        return myNumTaxa;
    }

    public int numSites() {
        return myNumSites;
    }

    public int numAlleles() {
        return myValues.size();
    }

    Collection<Byte2D> byteStorage() {
        return myValues.values();
    }
}
