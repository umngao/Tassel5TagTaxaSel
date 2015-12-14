/*
 *  DistanceMatrixBuilder
 * 
 *  Created on Nov 20, 2015
 */
package net.maizegenetics.taxa.distance;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;

/**
 *
 * @author Terry Casstevens
 */
public class DistanceMatrixBuilder {

    public static final String MATRIX_TYPE = "Matrix_Type";
    public static final String CENTERED_IBS_SUMPK = "Centered_IBS.SumPk";
    
    public static final String IBS_DISTANCE_MATRIX_TYPE = "IBS_Distance_Matrix";
    public static final String IBS_DISTANCE_MATRIX_NUM_ALLELES = "IBS_Distance_Matrix.NumAlleles";
    public static final String IBS_DISTANCE_MATRIX_AVE_TOTAL_SITES = "IBS_Distance_Matrix.AverageTotalSites";
    public static final String IBS_DISTANCE_MATRIX_TRUE_IBS = "IBS_Distance_Matrix.TrueIBS";

    private final int myNumTaxa;
    private final TaxaList myTaxa;
    private final double[][] myMatrix;
    private GeneralAnnotation myAnnotation = null;
    private final TaxaListBuilder myTaxaBuilder;

    private DistanceMatrixBuilder(int numTaxa, TaxaList taxa) {
        myTaxa = taxa;
        myNumTaxa = numTaxa;
        myMatrix = new double[myNumTaxa][myNumTaxa];
        if (myTaxa == null) {
            myTaxaBuilder = new TaxaListBuilder();
        } else {
            myTaxaBuilder = null;
        }
    }

    public static DistanceMatrixBuilder getInstance(TaxaList taxa) {
        return new DistanceMatrixBuilder(taxa.numberOfTaxa(), taxa);
    }

    public static DistanceMatrixBuilder getInstance(int numTaxa) {
        return new DistanceMatrixBuilder(numTaxa, null);
    }

    public void set(int x, int y, double value) {
        myMatrix[x][y] = myMatrix[y][x] = value;
    }

    public void addTaxon(Taxon taxon) {
        if (myTaxaBuilder == null) {
            throw new IllegalStateException("DistanceMatrixBuilder: addTaxon: this builder was given Taxa List at creation.");
        }
        myTaxaBuilder.add(taxon);
    }

    public DistanceMatrixBuilder annotation(GeneralAnnotation annotation) {
        myAnnotation = annotation;
        return this;
    }

    public DistanceMatrix build() {
        if (myTaxaBuilder == null) {
            return new DistanceMatrix(myMatrix, myTaxa, myAnnotation);
        } else {
            return new DistanceMatrix(myMatrix, myTaxaBuilder.build(), myAnnotation);
        }
    }

}
