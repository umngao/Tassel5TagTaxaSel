/*
 *  Byte3D
 */
package net.maizegenetics.dna.snp.byte3d;

/**
 * @author Terry Casstevens
 */
public interface Byte3D {

    public byte valueForAllele(int taxon, int site, int allele);

    public byte[][] valuesForAllSites(int taxon);

    public int numTaxa();

    public int numSites();
    
    public int numAlleles();

}
