/*
 *  Byte2D
 */
package net.maizegenetics.dna.snp.byte2d;

import net.maizegenetics.dna.snp.GenotypeTable;

/**
 * @author Terry Casstevens
 */
public interface Byte2D {

    public byte valueForAllele(int taxon, int site);

    public byte[] valuesForAllSites(int taxon);

    public int numTaxa();

    public int numSites();

    public GenotypeTable.SITE_SCORE_TYPE siteScoreType();

}
