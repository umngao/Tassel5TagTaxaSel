/*
 *  Dosage
 */
package net.maizegenetics.dna.snp.score;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class Dosage implements SiteScore {

    private final Byte2D myStorage;
    private final int myNumTaxa;
    private final int myNumSites;

    Dosage(Byte2D value) {
        myStorage = value;
        myNumTaxa = myStorage.numTaxa();
        myNumSites = myStorage.numSites();
    }

    Dosage(int numTaxa, int numSites) {
        myStorage = null;
        myNumTaxa = numTaxa;
        myNumSites = numSites;
    }

    public byte value(int taxon, int site) {
        return myStorage.valueForAllele(taxon, site);
    }

    Byte2D byteStorage() {
        return myStorage;
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        return new HashSet<>(Arrays.asList(SITE_SCORE_TYPE.Dosage));
    }

    @Override
    public int numTaxa() {
        return myNumTaxa;
    }

    @Override
    public int numSites() {
        return myNumSites;
    }
}
