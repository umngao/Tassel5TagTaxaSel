/*
 *  SiteScore
 */
package net.maizegenetics.dna.snp.score;

import java.util.Set;

/**
 *
 * @author Terry Casstevens
 */
public interface SiteScore {

    public static enum SITE_SCORE_TYPE {

        None, QualityScore, ReferenceProbablity, Dosage,
        DepthA, DepthC, DepthG, DepthT, DepthGap, DepthInsertion,
        ProbA, ProbC, ProbG, ProbT, ProbGap, ProbInsertion
    };

    /**
     * Return the site scores types.
     *
     * @return site score types.
     */
    public Set<SITE_SCORE_TYPE> siteScoreTypes();

    public int numTaxa();

    public int numSites();

//    public int numAlleles() {
//        return myValues.size();
//    }
//    protected Byte2D byteStorage(SITE_SCORE_TYPE type) {
//        return myValues.get(type);
//    }
//
//    Collection<Byte2D> byteStorage() {
//        return myValues.values();
//    }
}
