/*
 *  MemoryByte3D
 */
package net.maizegenetics.dna.snp.byte3d;

import net.maizegenetics.util.SuperByteMatrix;

/**
 * @author Terry Casstevens
 */
public class MemoryByte3D extends AbstractByte3D {

    private final SuperByteMatrix[] myValues;

    MemoryByte3D(SuperByteMatrix[] values, int numTaxa, int numSites) {
        super(values[0].getNumColumns(), numTaxa, numSites);
        myValues = values;
    }

    @Override
    public byte valueForAllele(int taxon, int site, int allele) {
        return myValues[taxon].get(site, allele);
    }
}
