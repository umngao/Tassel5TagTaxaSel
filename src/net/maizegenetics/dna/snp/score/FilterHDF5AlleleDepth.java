/*
 *  FilterHDF5AlleleDepth
 * 
 *  Created on May 2, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.FilterGenotypeTable;

/**
 *
 * @author Terry Casstevens
 */
public class FilterHDF5AlleleDepth extends AlleleDepth {

    private final HDF5AlleleDepth myBase;
    private final FilterGenotypeTable myFilter;

    FilterHDF5AlleleDepth(HDF5AlleleDepth alleleDepth, FilterGenotypeTable filterGenotypeTable) {
        super(filterGenotypeTable.numberOfTaxa(), filterGenotypeTable.numberOfSites());
        myBase = alleleDepth;
        myFilter = filterGenotypeTable;
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
        return myBase.valuesByte(myFilter.translateTaxon(taxon), myFilter.translateSite(site));
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
        if (!myFilter.isSiteFilterByRange() && !myFilter.isSiteFilter()) {
            return myBase.valuesForTaxonByte(myFilter.translateTaxon(taxon));
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
        return myBase.valueByte(myFilter.translateTaxon(taxon), myFilter.translateSite(site), scoreType);
    }

}
