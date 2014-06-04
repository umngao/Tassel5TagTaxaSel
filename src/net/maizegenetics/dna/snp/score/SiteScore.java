/*
 *  SiteScore
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.byte3d.Byte3D;

/**
 * In memory site scores class.
 *
 * @author Terry Casstevens
 */
public class SiteScore {

    private final GenotypeTable.SITE_SCORE_TYPE myScoreType;
    private final Byte3D myValues;

    SiteScore(Byte3D values, GenotypeTable.SITE_SCORE_TYPE scoreType) {
        myScoreType = scoreType;
        myValues = values;
    }

    public float siteScore(int taxon, int site, int allele) {
        return SiteScoreUtil.byteToFloatPercentage(myValues.valueForAllele(taxon, site, allele));
    }

    /**
     * Returns the site score of the given sequence and site. Default is allele
     * 0. Can be used with scores that only have one value per taxon and site.
     *
     * @param taxon taxon index
     * @param site site
     *
     * @return site score.
     */
    public float siteScore(int taxon, int site) {
        return siteScore(taxon, site, 0);
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
                result[t][s] = siteScore(t, s, 0);
            }
        }
        return result;
    }

    /**
     * Return what type of these site scores.
     *
     * @return site score type.
     */
    public GenotypeTable.SITE_SCORE_TYPE siteScoreType() {
        return myScoreType;
    }

    public int numTaxa() {
        return myValues.numTaxa();
    }

    public int numSites() {
        return myValues.numSites();
    }

    public int numAlleles() {
        return myValues.numAlleles();
    }
}
