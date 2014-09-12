/*
 *  AlleleDepth
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleDepth extends SiteScore {

    public static final SiteScore.SITE_SCORE_TYPE[] ALLELE_DEPTH_TYPES = new SiteScore.SITE_SCORE_TYPE[]{
        SiteScore.SITE_SCORE_TYPE.DepthA, SiteScore.SITE_SCORE_TYPE.DepthC,
        SiteScore.SITE_SCORE_TYPE.DepthG, SiteScore.SITE_SCORE_TYPE.DepthT,
        SiteScore.SITE_SCORE_TYPE.DepthGap, SiteScore.SITE_SCORE_TYPE.DepthInsertion};

    AlleleDepth(Byte2D[] values) {
        super(values);
    }

    public int value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        return AlleleDepthUtil.depthByteToInt(byteStorage(scoreType).valueForAllele(taxon, site));
    }

}
