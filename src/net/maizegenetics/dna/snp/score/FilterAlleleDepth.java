/*
 *  FilterAlleleDepth
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.Translate;

/**
 *
 * @author Terry Casstevens
 */
public class FilterAlleleDepth extends AlleleDepth {

    private final AlleleDepth myBase;
    private final Translate myTranslate;

    FilterAlleleDepth(AlleleDepth alleleDepth, Translate translate) {
        super(translate.numTaxa(), translate.numSites());
        myBase = alleleDepth;
        myTranslate = translate;
    }

    /**
     * Returns the depth values (byte representation) of all nucleotides at
     * given taxon and site. Depth values are stored in bytes and translated to
     * integer using AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     *
     * @return depths
     */
    @Override
    public byte[] valuesByte(int taxon, int site) {
        return myBase.valuesByte(myTranslate.taxon(taxon), myTranslate.site(site));
    }

    /**
     * Returns depth values (byte representation) of all nucleotides and sites
     * for given taxon. The first dimension of returned array is nucleotides
     * (ALLELE_DEPTH_TYPES) and second dimension is sites.
     *
     * @param taxon taxon
     *
     * @return depths
     */
    @Override
    public byte[][] valuesForTaxonByte(int taxon) {
        if (!myTranslate.hasSiteTranslations()) {
            return myBase.valuesForTaxonByte(myTranslate.taxon(taxon));
        } else {
            return super.valuesForTaxonByte(taxon);
        }
    }

    /**
     * Returns the depth (byte representation) of nucleotide (scoreType) at
     * given taxon and site. Depth values are stored in bytes and translated to
     * integer using AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     * @param scoreType nucleotide
     *
     * @return depth
     */
    @Override
    public byte valueByte(int taxon, int site, SiteScore.SITE_SCORE_TYPE scoreType) {
        return myBase.valueByte(myTranslate.taxon(taxon), myTranslate.site(site), scoreType);
    }

}

