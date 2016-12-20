/*
 *  TranslateIndex
 * 
 *  Created on Dec 14, 2016
 */
package net.maizegenetics.dna.snp;

/**
 * No translation to index.
 *
 * @author Terry Casstevens
 */
public class TranslateIndex {

    private final int myNumIndices;

    /**
     * Constructor
     *
     * @param numIndices number of sites
     */
    TranslateIndex(int numIndices) {
        myNumIndices = numIndices;
    }

    /**
     * Translates index to base index. This class has no translation.
     *
     * @param index index
     * @return translated base index
     */
    public int translate(int index) {
        return index;
    }

    /**
     * Translates base index to this index. This class has no translation.
     *
     * @param index index
     * @return translated index
     */
    public int reverseTranslateIndex(int index) {
        return index;
    }

    /**
     * Number of indices represented by this translation. Number of base indices
     * will be the same or larger.
     *
     * @return number of indices
     */
    public int numIndices() {
        return myNumIndices;
    }

    public boolean hasTranslations() {
        return false;
    }

    public int[] getTranslations() {
        int[] result = new int[numIndices()];
        for (int i = 0; i < numIndices(); i++) {
            result[i] = translate(i);
        }
        return result;
    }

}
