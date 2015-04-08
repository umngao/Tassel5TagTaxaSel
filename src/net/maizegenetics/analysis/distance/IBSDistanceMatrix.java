/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ProgressListener;

import java.util.concurrent.atomic.LongAdder;
import java.util.stream.IntStream;

import static net.maizegenetics.dna.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor;
import static net.maizegenetics.dna.WHICH_ALLELE.Minor2;

/**
 * This class calculates an identity by state matrix. It is scaled so only
 * non-missing comparison are used. It conducts bit level calculations of IBS
 * for genotypes. Only the two most common alleles are used in the distance
 * calculations.
 * <p>
 * Please note that when heterozygous genotypes are used, Het to Het distance is
 * 0.5 NOT 0.0. The default along the identity diagonal is 0 (isTrueIBS =
 * false), but changing isTrueIBS = true will calculate the identity.
 * <p>
 * The distance estimates become wildly inaccurate when too few sites are used
 * to calculate distance. The minSiteComp parameter can be used to control the
 * minimum number of sites used for a calculation. If there are insufficient
 * sites in the estimate, then Double.NaN is returned.
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class IBSDistanceMatrix extends DistanceMatrix {

    private ProgressListener myListener = null;
    private int numSeqs;
    private GenotypeTable theTBA = null;
    /**
     * Holds the average numbers of sites in the comparisons
     */
    private double avgTotalSites;
    private int minSitesComp = 0;
    private boolean isTrueIBS = false;

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     */
    public IBSDistanceMatrix(GenotypeTable theAlignment) {
        this(theAlignment, null);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     * @param listener Listener to track progress in calculations
     */
    public IBSDistanceMatrix(GenotypeTable theAlignment, ProgressListener listener) {
        this(theAlignment, 0, listener);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     * @param minSiteComp Minimum number of sites needed to estimate distance
     * @param listener Listener to track progress in calculations
     */
    public IBSDistanceMatrix(GenotypeTable theAlignment, int minSiteComp, ProgressListener listener) {
        this(theAlignment, minSiteComp, false, listener, true);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     * @param minSiteComp Minimum number of sites needed to estimate distance
     * @param trueIBS estimate diagonal distance based IBS (default = false,
     * i=i=0.0)
     * @param listener Listener to track progress in calculations
     * @param useThirdState
     */
    public IBSDistanceMatrix(GenotypeTable theAlignment, int minSiteComp, boolean trueIBS, ProgressListener listener, boolean useThirdState) {
        super();
        this.minSitesComp = minSiteComp;
        isTrueIBS = trueIBS;
        myListener = listener;
        numSeqs = theAlignment.numberOfTaxa();
        theTBA = theAlignment;
        //  this should have an option to only use the 2 or 3 most common alleles
        setIdGroup(theAlignment.taxa());
        computeHetBitDistances(useThirdState);
    }

    /**
     * This is a cleanest, fastest and most accurate way to calculate distance.
     */
    private void computeHetBitDistances(boolean useThirdState) {
        avgTotalSites = 0;
        LongAdder count = new LongAdder();
        //note this distance object is modified by a parallel stream, but each element is only touched once
        double[][] distance = new double[numSeqs][numSeqs];
        long numberOfTests=numSeqs*(numSeqs-1)/2;
        long time=System.currentTimeMillis();
        IntStream.range(0,numSeqs).parallel().forEach( i -> {
            long[] iMj = theTBA.allelePresenceForAllSites(i, Major).getBits();
            long[] iMn = theTBA.allelePresenceForAllSites(i, Minor).getBits();
            long[] iMn2 = null;
            if (useThirdState) {
                iMn2 = theTBA.allelePresenceForAllSites(i, Minor2).getBits();
            }
            for (int j = i; j < numSeqs; j++) {
                if (j == i && !isTrueIBS) {
                    distance[i][i] = 0;
                } else {
                    long[] jMj = theTBA.allelePresenceForAllSites(j, Major).getBits();
                    long[] jMn = theTBA.allelePresenceForAllSites(j, Minor).getBits();
                    double[] result;
                    if (useThirdState) {
                        long[] jMn2 = theTBA.allelePresenceForAllSites(j, Minor2).getBits();
                        result = computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesComp);
                    } else {
                        result = computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesComp);
                    }
                    distance[i][j] = distance[j][i] = result[0];
                    avgTotalSites += result[1];  //this assumes not hets
                    count.increment();
                }
            }
            fireProgress((int)((count.doubleValue() / (double)numberOfTests) * 50.0));
        });
        setDistances(distance);
        avgTotalSites /= (double) count.longValue();
        System.out.println("computeHetBitDistances time = " + (System.currentTimeMillis()-time)/1000 +" seconds");
    }

    /**
     * Compute distance for a pair of taxa.
     *
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(GenotypeTable theTBA, int taxon1, int taxon2) {
        return computeHetBitDistances(theTBA, taxon1, taxon2, 0);
    }

    /**
     * Compute distance for a pair of taxa.
     *
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     *
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(GenotypeTable theTBA, int taxon1, int taxon2, int minSitesCompared) {
        long[] iMj = theTBA.allelePresenceForAllSites(taxon1, Major).getBits();
        long[] iMn = theTBA.allelePresenceForAllSites(taxon1, Minor).getBits();
        long[] iMn2 = theTBA.allelePresenceForAllSites(taxon1, Minor2).getBits();
        long[] jMj = theTBA.allelePresenceForAllSites(taxon2, Major).getBits();
        long[] jMn = theTBA.allelePresenceForAllSites(taxon2, Minor).getBits();
        long[] jMn2 = theTBA.allelePresenceForAllSites(taxon2, Minor2).getBits();
        return computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesCompared, 0, iMj.length - 1);
    }

    /**
     * Compute distance for a pair of taxa. Optimized for calculations sites
     * within a certain range of underlying word (64 sites chunks) in the TBit
     * array
     *
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @param firstWord starting word for calculating distance
     * site=(firstWord*64)
     * @param lastWord ending word for calculating distance inclusive
     * site=(lastWord*64+63)
     * @param maskBadSet Optional mask for sites (those set to 1 are kept)
     *
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(GenotypeTable theTBA, int taxon1, int taxon2,
            int minSitesCompared, int firstWord, int lastWord, BitSet maskBadSet) {
        long[] iMj = theTBA.allelePresenceForAllSites(taxon1, Major).getBits();
        long[] iMn = theTBA.allelePresenceForAllSites(taxon1, Minor).getBits();
        long[] iMn2 = theTBA.allelePresenceForAllSites(taxon1, Minor2).getBits();
        if (maskBadSet != null) {
            long[] maskBad = maskBadSet.getBits();
            for (int i = 0; i < iMj.length; i++) {
                iMj[i] = iMj[i] & maskBad[i];
            }
            for (int i = 0; i < iMn.length; i++) {
                iMn[i] = iMn[i] & maskBad[i];
            }
        }
        long[] jMj = theTBA.allelePresenceForAllSites(taxon2, Major).getBits();
        long[] jMn = theTBA.allelePresenceForAllSites(taxon2, Minor).getBits();
        long[] jMn2 = theTBA.allelePresenceForAllSites(taxon2, Minor2).getBits();
        return computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesCompared, firstWord, lastWord);
    }

    /**
     * Calculation of distance using the bit vector of major and minor alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn, int minSitesCompared) {
        return computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesCompared, 0, iMj.length - 1);
    }

    /**
     * Calculation of distance using the bit vector of the first three alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistancesThirdState(long[] iMj, long[] iMn, long[] iMn2, long[] jMj, long[] jMn, long[] jMn2,
            int minSitesCompared) {
        return computeHetBitDistancesThirdState(iMj, iMn, iMn2, jMj, jMn, jMn2, minSitesCompared, 0, iMj.length - 1);
    }

    /**
     * Calculation of distance using the bit vector of major and minor alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @param firstWord first world for calculating distance
     * @param lastWord last word for calculating distance inclusive
     * site=(endWord*64+63)
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn,
            int minSitesCompared, int firstWord, int lastWord) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = firstWord; x <= lastWord; x++) {
            long same = (iMj[x] & jMj[x]) | (iMn[x] & jMn[x]);
            long diff = (iMj[x] & jMn[x]) | (iMn[x] & jMj[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        int sites = sameCnt + diffCnt - hetCnt;
        double identity = ((double) (sameCnt) - 0.5 * hetCnt) / (double) (sites);
        double dist = 1 - identity;
        if (sites > minSitesCompared) {
            return new double[]{dist, (double) sites};
        } else {
            return new double[]{Double.NaN, (double) sites};
        }
    }

    /**
     * Calculation of distance using the bit vectors of the first three alleles.
     *
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param iMn2 Vector of second minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param jMn2 Vector of second minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate
     * distance
     * @param firstWord first world for calculating distance
     * @param lastWord last word for calculating distance inclusive
     * site=(endWord*64+63)
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistancesThirdState(long[] iMj, long[] iMn, long[] iMn2, long[] jMj, long[] jMn, long[] jMn2,
            int minSitesCompared, int firstWord, int lastWord) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = firstWord; x <= lastWord; x++) {
            long same = (iMj[x] & jMj[x]) | (iMn[x] & jMn[x]) | (iMn2[x] & jMn2[x]);
            long diff = (iMj[x] & jMn[x]) | (iMn[x] & jMj[x]) | (iMj[x] & jMn2[x])
                    | (iMn2[x] & jMj[x]) | (iMn[x] & jMn2[x]) | (iMn2[x] & jMn[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        int sites = sameCnt + diffCnt - hetCnt;
        double identity = ((double) (sameCnt) - 0.5 * hetCnt) / (double) (sites);
        double dist = 1 - identity;
        if (sites > minSitesCompared) {
            return new double[]{dist, (double) sites};
        } else {
            return new double[]{Double.NaN, (double) sites};
        }
    }

    /*
     * Average number of sites used in calculating the distance matrix
     */
    public double getAverageTotalSites() {
        return avgTotalSites;
    }

    public String toString(int d) {
        double[][] distance = this.getDistances();
        /*Return a string representation of this matrix with 'd'
         displayed digits*/
        String newln = System.getProperty("line.separator");
        String outPut = new String();
        int i, j;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(5);
        for (i = 0; i < distance.length; i++) {
            for (j = 0; j < distance[i].length; j++) {

                //Numeric x = new Numeric(distance[i][j]);
                String num = nf.format(d);
                //num = x.toString(d);
                //ed change that screws up formatting
                //num=""+this.element[i][j];
                outPut = outPut + num + (char) 9;
            }
            outPut = outPut + newln;
        }
        return outPut;
    }

    @Override
    public String toString() {
        return this.toString(6);
    }

    protected void fireProgress(int percent) {
        if (myListener != null) {
            myListener.progress(percent, null);
        }

    }

    /*
     * Returns whether true IBS is calculated for the diagonal 
     */
    public boolean isTrueIBS() {
        return isTrueIBS;
    }
}
