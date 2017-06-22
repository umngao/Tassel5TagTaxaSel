/*
 *  MaskGenotypeCallTable
 * 
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.MaskMatrix;

/**
 * @author Terry Casstevens
 */
public class MaskGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable myBase;
    private final MaskMatrix myMask;

    public MaskGenotypeCallTable(GenotypeCallTable base, MaskMatrix mask) {
        super(base.numberOfTaxa(), base.numberOfSites(), base.isPhased(), base.alleleDefinitions());
        myBase = base;
        myMask = mask;
    }

    @Override
    public byte genotype(int taxon, int site) {
        if (myMask.get(taxon, site)) {
            return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        } else {
            return myBase.genotype(taxon, site);
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myBase.diploidAsString(site, genotype(taxon, site));
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(genotypeAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myBase.diploidAsString(site, value);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBase.transposeData(siteInnerLoop);
    }

    @Override
    public boolean isSiteOptimized() {
        return myBase.isSiteOptimized();
    }

}
