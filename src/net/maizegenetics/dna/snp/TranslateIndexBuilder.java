/*
 *  TranslateIndexBuilder
 * 
 *  Created on Dec 14, 2016
 */
package net.maizegenetics.dna.snp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateIndexBuilder {

    private final TranslateIndex myBase;
    private final int myNumBaseIndices;
    private boolean myHasNegativeIndices = false;
    private final List<Integer> myIndicesToKeep = new ArrayList<>();

    private TranslateIndexBuilder(int numBaseIndices) {
        myNumBaseIndices = numBaseIndices;
        myBase = null;
    }

    private TranslateIndexBuilder(TranslateIndex base) {
        myNumBaseIndices = base.numIndices();
        myBase = base;
    }

    public static TranslateIndexBuilder getInstance(int numBaseIndices) {
        return new TranslateIndexBuilder(numBaseIndices);
    }

    public static TranslateIndexBuilder getInstance(TranslateIndex base) {
        return new TranslateIndexBuilder(base);
    }

    public static TranslateIndexBuilder getInstance(int numBaseIndices, TranslateIndex base) {
        if (base == null) {
            return new TranslateIndexBuilder(numBaseIndices);
        } else {
            if (numBaseIndices != base.numIndices()) {
                throw new IllegalArgumentException("TranslateTaxaBuilder: getInstance: numBaseIndices: " + numBaseIndices + " should equal base: " + base.numIndices());
            }
            return new TranslateIndexBuilder(base);
        }
    }

    /**
     * Returns no translation instance.
     *
     * @param numIndices number of indices
     * @return no translation instance
     */
    public static TranslateIndex noTranslation(int numIndices) {
        return new TranslateIndex(numIndices);
    }

    public static TranslateIndex orderedTranslation(int[] indexRedirect, TranslateIndex base) {

        int numIndices = indexRedirect.length;

        if (base == null) {
            int[] result = new int[numIndices];
            for (int i = 0; i < numIndices; i++) {
                if (indexRedirect[i] < -1) {
                    throw new IllegalArgumentException("TranslateIndexBuilder: orderedTranslation: ordered translation can't contain negative indices.");
                }
                result[i] = indexRedirect[i];
            }
            Arrays.sort(result);
            return new TranslateIndexRedirect(result);
        } else {
            int[] result = new int[numIndices];
            for (int i = 0; i < numIndices; i++) {
                if (indexRedirect[i] < -1) {
                    throw new IllegalArgumentException("TranslateIndexBuilder: orderedTranslation: ordered translation can't contain negative indices.");
                }
                result[i] = base.translate(indexRedirect[i]);
            }
            Arrays.sort(result);
            return new TranslateIndexRedirect(result);
        }

    }

    public static TranslateIndex unorderedTranslation(int[] indicesNewOrder, TranslateIndex base) {

        int numIndices = indicesNewOrder.length;

        if (base == null) {
            int[] result = new int[numIndices];
            for (int i = 0; i < numIndices; i++) {
                if (indicesNewOrder[i] == -1) {
                    result[i] = -1;
                } else if (indicesNewOrder[i] < 0) {
                    throw new IllegalArgumentException("TranslateIndexBuilder: unorderedTranslation: only negative number allowed is -1.");
                } else {
                    result[i] = indicesNewOrder[i];
                }
            }
            return new TranslateIndexRedirectUnordered(result);
        } else {
            if (numIndices != base.numIndices()) {
                throw new IllegalStateException("TranslateIndexBuilder: unorderedTranslation: number of newly ordered indices: " + numIndices + " should equal base num indices: " + base.numIndices());
            }
            int[] result = new int[numIndices];
            for (int i = 0; i < numIndices; i++) {
                if (indicesNewOrder[i] == -1) {
                    result[i] = -1;
                } else if (indicesNewOrder[i] < 0) {
                    throw new IllegalArgumentException("TranslateIndexBuilder: unorderedTranslation: only negative number allowed is -1.");
                } else {
                    result[i] = base.translate(indicesNewOrder[i]);
                }
            }
            return new TranslateIndexRedirectUnordered(result);
        }

    }

    /**
     * Keeps a range of indices from start (inclusive) to end (inclusive)
     *
     * @param start start index
     * @param end end index
     * @param base base translation
     *
     * @return new translation
     */
    public static TranslateIndex range(int start, int end, TranslateIndex base) {

        if (start < 0 || start > end) {
            throw new IllegalArgumentException();
        }

        if (base == null) {
            return new TranslateIndexRange(start, end);
        } else {
            return new TranslateIndexBuilder(base).keepIndices(start, end).build();
        }

    }

    /**
     * Keep specified index.
     *
     * @param index index
     *
     * @return this builder
     */
    public TranslateIndexBuilder keepIndex(int index) {
        if (index == -1) {
            myHasNegativeIndices = true;
        } else if (index < 0 || index >= myNumBaseIndices) {
            throw new IllegalArgumentException();
        }
        myIndicesToKeep.add(index);
        return this;
    }

    /**
     * Keeps a range of indices from start (inclusive) to end (inclusive)
     *
     * @param start lower index
     * @param end last index to keep
     *
     * @return this builder
     */
    public TranslateIndexBuilder keepIndices(int start, int end) {
        if (start < 0 || start > end || end >= myNumBaseIndices) {
            throw new IllegalArgumentException();
        }
        for (int i = start; i <= end; i++) {
            myIndicesToKeep.add(i);
        }
        return this;
    }

    public TranslateIndexBuilder keepIndices(int[] indices) {
        for (int current : indices) {
            keepIndex(current);
        }
        return this;
    }

    public int numIndices() {
        return myIndicesToKeep.size();
    }

    public TranslateIndex build() {

        int numIndicesToKeep = myIndicesToKeep.size();

        if (numIndicesToKeep == 0) {
            throw new IllegalStateException("TranslateIndexBuilder: build: no indices to keep.");
        } else if (numIndicesToKeep == myNumBaseIndices && !myHasNegativeIndices) {
            if (myBase != null) {
                return myBase;
            } else {
                return new TranslateIndex(myNumBaseIndices);
            }
        }

        int[] indexRedirect = new int[numIndicesToKeep];
        if (myBase == null) {
            for (int i = 0; i < numIndicesToKeep; i++) {
                indexRedirect[i] = myIndicesToKeep.get(i);
            }
        } else {
            for (int i = 0; i < numIndicesToKeep; i++) {
                indexRedirect[i] = myBase.translate(myIndicesToKeep.get(i));
            }
        }

        if (myHasNegativeIndices) {
            return new TranslateIndexRedirectUnordered(indexRedirect);
        } else {
            Arrays.sort(indexRedirect);
            return new TranslateIndexRedirect(indexRedirect);
        }

    }

}
