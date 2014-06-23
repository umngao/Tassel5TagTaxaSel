/*
 *  AlleleDepth
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleDepth extends SiteScore {

    AlleleDepth(Byte2D[] values) {
        super(values);
    }

    public int value(int taxon, int site, GenotypeTable.SITE_SCORE_TYPE scoreType) {
        return AlleleDepthUtil.depthByteToInt(byteStorage(scoreType).valueForAllele(taxon, site));
    }

}
