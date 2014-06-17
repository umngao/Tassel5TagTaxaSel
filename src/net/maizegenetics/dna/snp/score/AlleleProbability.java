/*
 *  AlleleProbability
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleProbability extends SiteScore {

    AlleleProbability(Byte2D[] values) {
        super(values);
    }

    public float value(int taxon, int site, GenotypeTable.SITE_SCORE_TYPE scoreType) {
        return SiteScoreUtil.byteToFloatPercentage(byteStorage(scoreType).valueForAllele(taxon, site));
    }

}
