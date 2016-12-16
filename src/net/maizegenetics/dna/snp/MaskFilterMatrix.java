/*
 *  MaskFilterMatrix
 * 
 *  Created on Dec 13, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Terry Casstevens
 */
public class MaskFilterMatrix implements MaskMatrix {

    private final MaskMatrix myBase;
    private final Translate myTranslate;

    MaskFilterMatrix(MaskMatrix base, Translate translate) {
        if (base instanceof MaskFilterMatrix) {
            throw new IllegalArgumentException();
        }
        myBase = base;
        myTranslate = translate;
    }

    @Override
    public boolean get(int taxon, int site) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return true;
        } else {
            return myBase.get((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF));
        }
    }

    @Override
    public boolean isTaxonMaskedHint(int taxon) {
        int newTaxon = myTranslate.taxon(taxon);
        if (newTaxon == -1) {
            return true;
        } else {
            return myBase.isTaxonMaskedHint(newTaxon);
        }
    }

    @Override
    public BitSet maskForTaxon(int taxon) {
        int newTaxon = myTranslate.taxon(taxon);
        if (newTaxon == -1) {
            BitSet result = new OpenBitSet(numSites());
            result.set(0, numSites());
            return result;
        } else {
            return myBase.maskForTaxon(newTaxon);
        }
    }

    @Override
    public boolean isSiteMaskedHint(int site) {
        int newSite = myTranslate.site(site);
        if (newSite == -1) {
            return true;
        } else {
            return myBase.isSiteMaskedHint(newSite);
        }
    }

    @Override
    public BitSet maskForSite(int site) {
        int newSite = myTranslate.site(site);
        if (newSite == -1) {
            BitSet result = new OpenBitSet(numTaxa());
            result.set(0, numTaxa());
            return result;
        } else {
            return myBase.maskForSite(newSite);
        }
    }

    @Override
    public int numTaxa() {
        return myTranslate.numTaxa();
    }

    @Override
    public int numSites() {
        return myTranslate.numSites();
    }

    @Override
    public boolean isSiteOptimized() {
        return myBase.isSiteOptimized();
    }

}
