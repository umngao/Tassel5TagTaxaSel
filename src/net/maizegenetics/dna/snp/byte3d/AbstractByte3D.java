/*
 *  AbstractByte3D
 */
package net.maizegenetics.dna.snp.byte3d;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractByte3D implements Byte3D {

    private final int myNumAlleles;
    private final int myNumTaxa;
    private final int myNumSites;

    public AbstractByte3D(int numAlleles, int numTaxa, int numSites) {
        myNumAlleles = numAlleles;
        myNumTaxa = numTaxa;
        myNumSites = numSites;
    }

    @Override
    public byte[][] valuesForAllSites(int taxon) {
        byte[][] result = new byte[myNumAlleles][myNumSites];
        for (int allele = 0; allele < myNumAlleles; allele++) {
            for (int site = 0; site < myNumSites; site++) {
                result[allele][site] = valueForAllele(taxon, site, allele);
            }
        }
        return result;
    }

    @Override
    public int numTaxa() {
        return myNumTaxa;
    }

    @Override
    public int numSites() {
        return myNumSites;
    }

    @Override
    public int numAlleles() {
        return myNumAlleles;
    }

}
