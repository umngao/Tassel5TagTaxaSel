/*
 *  MemoryByte2D
 */
package net.maizegenetics.dna.snp.byte2d;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.SuperByteMatrix;

/**
 * @author Terry Casstevens
 */
public class MemoryByte2D extends AbstractByte2D {

    private final SuperByteMatrix myValues;

    MemoryByte2D(GenotypeTable.SITE_SCORE_TYPE scoreType, SuperByteMatrix values) {
        super(scoreType, values.getNumRows(), values.getNumColumns());
        myValues = values;
    }

    @Override
    public byte valueForAllele(int taxon, int site) {
        return myValues.get(taxon, site);
    }
}
