/*
 *  DistanceMatrixWithCounts
 * 
 *  Created on Jan 14, 2016
 */
package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.GeneralAnnotation;

/**
 *
 * @author Terry Casstevens
 */
public class DistanceMatrixWithCounts extends DistanceMatrix {

    private final int[] myCounts;

    DistanceMatrixWithCounts(double[][] distance, TaxaList taxaList, GeneralAnnotation annotations, int[] counts) {
        super(distance, taxaList, annotations);
        myCounts = counts;
    }

    public int getCount(int x, int y) {

        if (x > y) {
            int temp = y;
            y = x;
            x = temp;
        }

        int index = y * (y + 1) / 2 + x;

        return myCounts[index];

    }

}
