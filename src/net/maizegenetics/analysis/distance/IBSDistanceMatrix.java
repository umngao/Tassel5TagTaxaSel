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

import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ProgressListener;

import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

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
    private final int numSeqs;
    private final GenotypeTable theTBA;
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
        //setupHetBitDistancesIncrementors();
        //computeHetBitDistances();
    }

    private void computeHetBitDistances() {
        avgTotalSites = 0;
        long time = System.currentTimeMillis();

        int[][] counters = stream(theTBA).collect(toCounters(numSeqs)).myCounters;

        double[][] distance = new double[numSeqs][numSeqs];
        long count = 0;
        for (int i = 0; i < numSeqs; i++) {
            int index = 0;
            for (int j = i; j < numSeqs; j++) {
                if (j == i && !isTrueIBS) {
                    distance[i][i] = 0;
                    index += 3;
                } else {
                    int sameCount = counters[i][index++];
                    int diffCount = counters[i][index++];
                    int hetCount = counters[i][index++];
                    long sites = sameCount + diffCount - hetCount;
                    double identity = ((double) (sameCount) - 0.5 * hetCount) / (double) (sites);
                    double dist = 1 - identity;

                    if (sites < minSitesComp) {
                        dist = Double.NaN;
                    }
                    distance[i][j] = distance[j][i] = dist;
                    avgTotalSites += sites;  //this assumes not hets
                    count++;
                }
            }
        }
        setDistances(distance);

        avgTotalSites /= (double) count;
        System.out.println("computeHetBitDistances time = " + (System.currentTimeMillis() - time) / 1000 + " seconds");
    }

    public static double[] computeHetDistances(byte[] first, byte[] second, int minSitesComp) {
        setupHetBitDistancesIncrementors();
        long counts = 0;
        for (int i = 0; i < first.length; i++) {
            int key1 = (first[i] & 0x70) >>> 1 | (first[i] & 0x7);
            int key2 = (second[i] & 0x70) >>> 1 | (second[i] & 0x7);
            counts += INCREMENT_FUNCTIONS[key1][key2];
        }

        int sameCount = (int) (counts & 0x1FFFFFl);
        int diffCount = (int) ((counts >>> 21) & 0x1FFFFFl);
        int hetCount = (int) ((counts >>> 42) & 0x1FFFFFl);
        long sites = sameCount + diffCount - hetCount;
        double identity = ((double) (sameCount) - 0.5 * hetCount) / (double) (sites);
        double dist = 1 - identity;

        if (sites < minSitesComp) {
            dist = Double.NaN;
        }

        return new double[]{dist, sites};
    }

    /**
     * This is a cleanest, fastest and most accurate way to calculate distance.
     */
    private void computeHetBitDistances(boolean useThirdState) {
        avgTotalSites = 0;
        //LongAdder count = new LongAdder();
        //note this distance object is modified by a parallel stream, but each element is only touched once
        double[][] distance = new double[numSeqs][numSeqs];
        long numberOfTests = numSeqs * (numSeqs - 1) / 2;
        long time = System.currentTimeMillis();
        IntStream.range(0, numSeqs).parallel().forEach(i -> {
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
                    //count.increment();
                }
            }
            fireProgress((int) ((double) (i + 1) / (double) numSeqs * 100.0));
        });
        setDistances(distance);
        //avgTotalSites /= (double) count.longValue();
        System.out.println("computeHetBitDistances time = " + (System.currentTimeMillis() - time) / 1000 + " seconds");
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

    private static final long TRUE_TRUE_LONG = 0x40000200001l;
    private static final long TRUE_FALSE_LONG = 0x1l;
    private static final long FALSE_TRUE_LONG = 0x200000l;
    private static final long FALSE_FALSE_LONG = 0x0l;
    private static long[][] INCREMENT_FUNCTIONS = null;

    public static void setupHetBitDistancesIncrementors() {

        if (INCREMENT_FUNCTIONS != null) {
            return;
        }

        INCREMENT_FUNCTIONS = new long[64][64];

        int[] possibleValues = new int[]{0, 1, 2, 3, 4, 5, 15};

        for (int a = 0; a < 7; a++) {
            for (int b = 0; b < 7; b++) {
                for (int c = 0; c < 7; c++) {
                    for (int d = 0; d < 7; d++) {

                        int key1 = (possibleValues[a] & 0x7) << 3 | (possibleValues[b] & 0x7);
                        int key2 = (possibleValues[c] & 0x7) << 3 | (possibleValues[d] & 0x7);

                        byte[] target = new byte[2];
                        target[0] = (byte) possibleValues[a];
                        target[1] = (byte) possibleValues[b];

                        byte[] match = new byte[2];
                        match[0] = (byte) possibleValues[c];
                        match[1] = (byte) possibleValues[d];

                        int[] result = new int[3];
                        result[0] = 0;
                        result[1] = 0;

                        if (target[0] != GenotypeTable.UNKNOWN_ALLELE) {
                            if ((target[0] == match[0]) || (target[0] == match[1])) {
                                result[0] = 1;
                            }
                            if ((match[0] != GenotypeTable.UNKNOWN_ALLELE) && (match[0] != target[0])) {
                                result[1] = 1;
                            } else if ((match[1] != GenotypeTable.UNKNOWN_ALLELE) && (match[1] != target[0])) {
                                result[1] = 1;
                            }
                        }
                        if ((result[0] == 1) && (result[1] == 1)) {
                            INCREMENT_FUNCTIONS[key1][key2] = TRUE_TRUE_LONG;
                            continue;
                        }
                        if (target[1] != GenotypeTable.UNKNOWN_ALLELE) {
                            if ((target[1] == match[0]) || (target[1] == match[1])) {
                                result[0] = 1;
                            }
                            if ((match[0] != GenotypeTable.UNKNOWN_ALLELE) && (match[0] != target[1])) {
                                result[1] = 1;
                            } else if ((match[1] != GenotypeTable.UNKNOWN_ALLELE) && (match[1] != target[1])) {
                                result[1] = 1;
                            }
                        }

                        if ((result[0] == 1) && (result[1] == 1)) {
                            INCREMENT_FUNCTIONS[key1][key2] = TRUE_TRUE_LONG;
                            continue;
                        }

                        if (result[0] == 1) {
                            INCREMENT_FUNCTIONS[key1][key2] = TRUE_FALSE_LONG;
                        } else if (result[1] == 1) {
                            INCREMENT_FUNCTIONS[key1][key2] = FALSE_TRUE_LONG;
                        } else {
                            INCREMENT_FUNCTIONS[key1][key2] = FALSE_FALSE_LONG;
                        }

                    }
                }
            }
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

    private static class Counters {

        private final int[][] myCounters;
        private final int myNumTaxa;

        public Counters(int numTaxa) {
            myNumTaxa = numTaxa;
            myCounters = new int[myNumTaxa][];
            for (int i = 0; i < myNumTaxa; i++) {
                myCounters[i] = new int[(myNumTaxa - i) * 3];
            }
        }

        public void add(byte[] values) {
            int taxaIndex = myNumTaxa - values.length;
            int[] currentCounter = myCounters[taxaIndex];
            int index = 0;
            for (int i = 0; i < values.length; i++) {
                currentCounter[index] += values[i] & 0x1;
                currentCounter[index + 1] += (values[i] >>> 1) & 0x1;
                currentCounter[index + 2] += (values[i] >>> 2) & 0x1;
                index += 3;
            }
        }

        public void add(long[] values) {
            int index = 0;
            for (int i = 0; i < myNumTaxa; i++) {
                for (int j = 0; j < myCounters[i].length; j += 3) {
                    myCounters[i][j] += (int) (values[index] & 0x1FFFFFl);
                    myCounters[i][j + 1] += (int) ((values[index] >>> 21) & 0x1FFFFFl);
                    myCounters[i][j + 2] += (int) ((values[index] >>> 42) & 0x1FFFFFl);
                    index++;
                }
            }
        }

        public void addAll(Counters counters) {
            int[][] other = counters.myCounters;
            for (int t = 0; t < myNumTaxa; t++) {
                for (int i = 0, n = myCounters[t].length; i < n; i++) {
                    myCounters[t][i] += other[t][i];
                }
            }
        }

    }

    private static Collector<long[], ?, Counters> toCounters(int numTaxa) {
        return new CountersCollector(numTaxa);
    }

    private static class CountersCollector implements Collector<long[], Counters, Counters> {

        private final int myNumTaxa;

        public CountersCollector(int numTaxa) {
            myNumTaxa = numTaxa;
        }

        @Override
        public Supplier<Counters> supplier() {
            return () -> new Counters(myNumTaxa);
        }

        @Override
        public BiConsumer<Counters, long[]> accumulator() {
            return Counters::add;
        }

        @Override
        public BinaryOperator<Counters> combiner() {
            return (left, right) -> {
                left.addAll(right);
                return left;
            };
        }

        @Override
        public Function<Counters, Counters> finisher() {
            return (Counters t) -> t;
        }

        @Override
        public Set<Collector.Characteristics> characteristics() {
            return Collections.unmodifiableSet(EnumSet.of(Collector.Characteristics.UNORDERED, Collector.Characteristics.IDENTITY_FINISH));
            //return Collections.EMPTY_SET;
        }

    }

    private int myNumSitesProcessed = 0;
    private static final int MAX_NUMBER_20_BITS = 0xFFFFF;

    public Stream<long[]> stream(GenotypeTable genotypes) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new ByteCounterSpliterator(genotypes, 0, genotypes.numberOfSites()), true);
    }

    class ByteCounterSpliterator implements Spliterator<long[]> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private byte[] myCachedSiteGenotype;
        private int myCachedSite;
        private final int myNumTaxa;
        private final int myNumSites;
        private int myFirstTaxa;

        ByteCounterSpliterator(GenotypeTable genotypes, int currentIndex, int fence) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myFirstTaxa = 0;
            myCurrentSite = currentIndex;
            myFence = fence;
            myCachedSiteGenotype = myGenotypes.genotypeAllTaxa(myCurrentSite);
            myCachedSite = myCurrentSite;
        }

        @Override
        public void forEachRemaining(Consumer<? super long[]> action) {
            int numSitesProcessed = myFence - myCurrentSite;
            for (; myCurrentSite < myFence; myCurrentSite += MAX_NUMBER_20_BITS) {
                int currentBlockFence = Math.min(myCurrentSite + MAX_NUMBER_20_BITS, myFence);
                long[] result = new long[myNumTaxa * (myNumTaxa + 1) / 2];
                for (; myCurrentSite < currentBlockFence; myCurrentSite++) {
                    if (myCurrentSite != myCachedSite) {
                        myCachedSiteGenotype = myGenotypes.genotypeAllTaxa(myCurrentSite);
                        myCachedSite = myCurrentSite;
                    }
                    int[] keys = new int[myNumTaxa];
                    for (int i = 0; i < myNumTaxa; i++) {
                        keys[i] = (myCachedSiteGenotype[i] & 0x70) >>> 1 | (myCachedSiteGenotype[i] & 0x7);
                    }
                    int index = 0;
                    for (; myFirstTaxa < myNumTaxa;) {
                        if (myCachedSiteGenotype[myFirstTaxa] != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                            long[] incrementFunctions = INCREMENT_FUNCTIONS[keys[myFirstTaxa]];
                            for (int secondTaxa = myFirstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                result[index++] += incrementFunctions[keys[secondTaxa]];
                            }
                        } else {
                            index += myNumTaxa - myFirstTaxa;
                        }
                        myFirstTaxa++;
                    }
                    myFirstTaxa = 0;
                }
                action.accept(result);
            }
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0));
        }

        @Override
        public boolean tryAdvance(Consumer<? super long[]> action) {
            if ((myCurrentSite < myFence) && (myFirstTaxa < myNumTaxa)) {
                int key1 = (myCachedSiteGenotype[myFirstTaxa] & 0x70) >>> 1 | (myCachedSiteGenotype[myFirstTaxa] & 0x7);
                long[] incrementFunctions = INCREMENT_FUNCTIONS[key1];

                long[] result = new long[myNumTaxa - myFirstTaxa];
                int index = 0;
                for (int secondTaxa = myFirstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                    int key2 = (myCachedSiteGenotype[secondTaxa] & 0x70) >>> 1 | (myCachedSiteGenotype[secondTaxa] & 0x7);
                    result[index++] = incrementFunctions[key2];
                }
                action.accept(result);

                myFirstTaxa++;
                if (myFirstTaxa == myNumTaxa) {
                    myCurrentSite++;
                    if (myCurrentSite < myFence) {
                        myCachedSiteGenotype = myGenotypes.genotypeAllTaxa(myCurrentSite);
                        myCachedSite = myCurrentSite;
                    } else {
                        myCachedSiteGenotype = null;
                        myCachedSite = -1;
                    }
                    myFirstTaxa = 0;
                }
                return true;
            } else {
                return false;
            }
        }

        @Override
        public Spliterator<long[]> trySplit() {
            int lo = myCurrentSite;
            int mid = (lo + myFence) >>> 1;
            if (lo < mid) {
                myCurrentSite = mid;
                return new ByteCounterSpliterator(myGenotypes, lo, mid);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            //return (long) (myFence - myCurrentSite) * (long) myNumTaxa;
            return (long) (myFence - myCurrentSite) / (long) MAX_NUMBER_20_BITS + 1l;
        }

        @Override
        public int characteristics() {
            return IMMUTABLE;
        }
    }

}
