/*
 *  TranslateSiteRange
 * 
 *  Created on May 7, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateSiteRange extends TranslateSite {

    private final int myRangeStart;
    private final int myRangeEnd;

    /**
     * Site Range Translation
     *
     * @param start start site (inclusive)
     * @param end end site (exclusive)
     */
    TranslateSiteRange(int start, int end) {
        super(end - start);
        myRangeStart = start;
        myRangeEnd = end;
    }

    @Override
    public int translateSite(int site) {
        return site + myRangeStart;
    }

    @Override
    public int reverseTranslateSite(int site) {
        return site - myRangeStart;
    }

}
