/*
 *  TranslateTaxaBuilder
 * 
 *  Created on May 27, 2016
 */
package net.maizegenetics.dna.snp;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Set;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateTaxaBuilder {
    
    private static final Logger myLogger = Logger.getLogger(TranslateTaxaBuilder.class);
    
    private final Set<Integer> myTaxaToKeep = new LinkedHashSet<>();
    private final TranslateTaxa myBase;
    private final int myNumBaseTaxa;
    private boolean myHasNegativeIndices = false;
    private boolean mySortAsOriginal = false;
    
    private TranslateTaxaBuilder(int numBaseTaxa) {
        myNumBaseTaxa = numBaseTaxa;
        myBase = null;
    }
    
    private TranslateTaxaBuilder(TranslateTaxa base) {
        myNumBaseTaxa = base.numTaxa();
        myBase = base;
    }
    
    public static TranslateTaxaBuilder getInstance(int numBaseTaxa) {
        return new TranslateTaxaBuilder(numBaseTaxa);
    }
    
    public static TranslateTaxaBuilder getInstance(TranslateTaxa base) {
        return new TranslateTaxaBuilder(base);
    }
    
    public static TranslateTaxaBuilder getInstance(int numBaseTaxa, TranslateTaxa base) {
        if (base == null) {
            return new TranslateTaxaBuilder(numBaseTaxa);
        } else {
            if (numBaseTaxa != base.numTaxa()) {
                throw new IllegalArgumentException("TranslateTaxaBuilder: getInstance: numBaseTaxa: " + numBaseTaxa + " should equal base: " + base.numTaxa());
            }
            return new TranslateTaxaBuilder(base);
        }
    }

    /**
     * Returns no translation instance.
     *
     * @param numTaxa number of taxa
     * @return no translation instance
     */
    public static TranslateTaxa noTranslation(int numTaxa) {
        return new TranslateTaxa(numTaxa);
    }
    
    public TranslateTaxaBuilder keepTaxon(int taxon) {
        if (taxon == -1) {
            myHasNegativeIndices = true;
        } else if (taxon >= myNumBaseTaxa) {
            throw new IllegalArgumentException("TranslateTaxaBuilder: keepTaxon: taxon: " + taxon + " must be less than: " + myNumBaseTaxa);
        }
        myTaxaToKeep.add(taxon);
        return this;
    }
    
    public TranslateTaxaBuilder sortAsOriginal() {
        mySortAsOriginal = true;
        return this;
    }
    
    public int numTaxa() {
        return myTaxaToKeep.size();
    }
    
    public TranslateTaxa build() {
        return build(false).x;
    }
    
    public Tuple<TranslateTaxa, int[]> build(boolean getOrigRedirect) {
        
        int numTaxaToKeep = myTaxaToKeep.size();
        
        if (numTaxaToKeep == 0) {
            throw new IllegalStateException("TranslateTaxaBuilder: build: no taxa to keep.");
        } else if (numTaxaToKeep == myNumBaseTaxa && !myHasNegativeIndices) {
            if (mySortAsOriginal) {
                if (myBase != null) {
                    return new Tuple<>(myBase, null);
                } else {
                    return new Tuple<>(new TranslateTaxa(myNumBaseTaxa), null);
                }
            } else {
                boolean sorted = true;
                int i = 0;
                for (Integer current : myTaxaToKeep) {
                    if (current != i) {
                        sorted = false;
                        break;
                    }
                    i++;
                }
                if (sorted) {
                    if (myBase != null) {
                        return new Tuple<>(myBase, null);
                    } else {
                        return new Tuple<>(new TranslateTaxa(myNumBaseTaxa), null);
                    }
                }
            }
        }
        
        int[] taxaRedirect = new int[numTaxaToKeep];
        int count = 0;
        for (Integer current : myTaxaToKeep) {
            taxaRedirect[count++] = current;
        }
        
        if (mySortAsOriginal) {
            if (myHasNegativeIndices) {
                myLogger.warn("build: sorting with negative indices.");
            }
            Arrays.sort(taxaRedirect);
        }
        
        int[] origRedirect = null;
        if (getOrigRedirect) {
            origRedirect = new int[numTaxaToKeep];
            System.arraycopy(taxaRedirect, 0, origRedirect, 0, numTaxaToKeep);
        }
        
        if (myBase != null) {
            for (int i = 0; i < numTaxaToKeep; i++) {
                taxaRedirect[i] = myBase.translate(taxaRedirect[i]);
            }
        }
        
        return new Tuple<>(new TranslateTaxaRedirect(taxaRedirect), origRedirect);
        
    }
    
}
