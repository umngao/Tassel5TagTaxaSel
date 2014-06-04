/*
 *  FilterByte3D
 */
package net.maizegenetics.dna.snp.byte3d;

import net.maizegenetics.dna.snp.FilterGenotypeTable;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class FilterByte3D extends AbstractByte3D {

    private final FilterGenotypeTable myFilterGenotypeTable;
    private final Byte3D myBaseByte3D;

    public FilterByte3D(Byte3D baseByte3D, FilterGenotypeTable filterGenotypeTable) {
        super(filterGenotypeTable.maxNumAlleles(), filterGenotypeTable.numberOfTaxa(), filterGenotypeTable.numberOfSites());
        myBaseByte3D = baseByte3D;
        myFilterGenotypeTable = filterGenotypeTable;
    }

    @Override
    public byte valueForAllele(int taxon, int site, int allele) {
        return myBaseByte3D.valueForAllele(myFilterGenotypeTable.translateTaxon(taxon),
                myFilterGenotypeTable.translateSite(site), allele);
    }
}
