/*
 *  GCTADistanceMatrix
 * 
 *  Created on May 31, 2015
 */
package net.maizegenetics.analysis.distance;

import java.util.Optional;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author Terry Casstevens
 */
public class GCTADistanceMatrix {

    private GCTADistanceMatrix() {
        // utility
    }

    /**
     * Compute GCTA kinship for all pairs of taxa. Missing sites are ignored.
     * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/pdf/main.pdf
     * Equation-3
     *
     * @param genotype Genotype Table used to compute kinship
     */
    public static DistanceMatrix getInstance(GenotypeTable genotype) {
        return getInstance(genotype, null);
    }

    public static DistanceMatrix getInstance(GenotypeTable genotype, ProgressListener listener) {
        return computeGCTADistances(genotype, listener);
    }

    private static DistanceMatrix computeGCTADistances(GenotypeTable genotype, ProgressListener listener) {

        int numSeqs = genotype.numberOfTaxa();
        long time = System.currentTimeMillis();

        Optional<CountersDistances> optional = stream(genotype, listener).reduce((CountersDistances t, CountersDistances u) -> {
            t.addAll(u);
            return t;
        });
        if (!optional.isPresent()) {
            return null;
        }

        CountersDistances counters = optional.get();
        int[] counts = counters.myCounters;
        float[] distances = counters.myDistances;

        double[][] result = new double[numSeqs][numSeqs];

        int index = 0;
        for (int t = 0; t < numSeqs; t++) {
            for (int i = 0, n = numSeqs - t; i < n; i++) {
                result[t][t + i] = result[t + i][t] = distances[index] / (double) counts[index];
                index++;
            }
        }

        System.out.println("computeGCTADistances time = " + (System.currentTimeMillis() - time) / 1000 + " seconds");
        return new DistanceMatrix(result, genotype.taxa());

    }

    protected static void fireProgress(int percent, ProgressListener listener) {
        if (listener != null) {
            listener.progress(percent, null);
        }

    }

    private static class CountersDistances {

        private final int[] myCounters;
        private final float[] myDistances;
        private final int myNumTaxa;

        public CountersDistances(int numTaxa) {
            myNumTaxa = numTaxa;
            myCounters = new int[myNumTaxa * (myNumTaxa + 1) / 2];
            myDistances = new float[myNumTaxa * (myNumTaxa + 1) / 2];
        }

        public void addAll(CountersDistances counters) {
            float[] otherDistances = counters.myDistances;
            for (int t = 0, n = myCounters.length; t < n; t++) {
                myDistances[t] += otherDistances[t];
            }
            otherDistances = null;
            int[] otherCounters = counters.myCounters;
            for (int t = 0, n = myCounters.length; t < n; t++) {
                myCounters[t] += otherCounters[t];
            }
        }

    }

    private static final byte[] PRECALCULATED_COUNTS = new byte[512];

    static {
        for (int major = 0; major < 8; major++) {
            for (int a = 0; a < 8; a++) {
                for (int b = 0; b < 8; b++) {
                    int temp = (major << 6) | (a << 3) | b;
                    if ((major == 7) | ((a == 7) && (b == 7))) {
                        PRECALCULATED_COUNTS[temp] = 3;
                    } else {
                        if (a == major) {
                            PRECALCULATED_COUNTS[temp]++;
                        }
                        if (b == major) {
                            PRECALCULATED_COUNTS[temp]++;
                        }
                    }
                }
            }
        }
    }

    private static final byte[] INCREMENT = new byte[4096];

    static {
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                int temp = a << 10 | b << 8;
                for (int c = 0; c < 4; c++) {
                    for (int d = 0; d < 4; d++) {
                        int temp2 = c << 6 | d << 4;
                        for (int e = 0; e < 4; e++) {
                            for (int f = 0; f < 4; f++) {
                                int incrementIndex = temp | temp2 | e << 2 | f;
                                if ((a != 3) && (b != 3)) {
                                    INCREMENT[incrementIndex]++;
                                }
                                if ((c != 3) && (d != 3)) {
                                    INCREMENT[incrementIndex]++;
                                }
                                if ((e != 3) && (f != 3)) {
                                    INCREMENT[incrementIndex]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private static int myNumSitesProcessed = 0;

    private static Stream<CountersDistances> stream(GenotypeTable genotypes, ProgressListener listener) {
        myNumSitesProcessed = 0;
        return StreamSupport.stream(new GCTASiteSpliterator(genotypes, 0, genotypes.numberOfSites(), listener), true);
    }

    static class GCTASiteSpliterator implements Spliterator<CountersDistances> {

        private int myCurrentSite;
        private final int myFence;
        private final GenotypeTable myGenotypes;
        private final int myNumTaxa;
        private final int myNumSites;
        private final ProgressListener myProgressListener;

        GCTASiteSpliterator(GenotypeTable genotypes, int currentIndex, int fence, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myCurrentSite = currentIndex;
            myFence = fence;
            myProgressListener = listener;
        }

        @Override
        public void forEachRemaining(Consumer<? super CountersDistances> action) {

            int numSitesProcessed = myFence - myCurrentSite;
            CountersDistances result = new CountersDistances(myNumTaxa);
            int[] counts = result.myCounters;
            float[] distances = result.myDistances;;

            for (; myCurrentSite < myFence; myCurrentSite += 3) {

                int numSitesPerBlock = Math.min(3, myFence - myCurrentSite);

                short[] majorCount = new short[myNumTaxa];
                float[] answer = new float[4096];

                byte major = myGenotypes.majorAllele(myCurrentSite);
                float majorFreq = (float) myGenotypes.majorAlleleFrequency(myCurrentSite);
                float majorFreqTimes2 = majorFreq * 2.0f;
                float denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                    float zeroTerm = 0.0f - majorFreqTimes2;
                    float oneTerm = 1.0f - majorFreqTimes2;
                    float twoTerm = 2.0f - majorFreqTimes2;

                    float zeroZero = zeroTerm * zeroTerm / denominatorTerm;
                    float zeroOne = zeroTerm * oneTerm / denominatorTerm;
                    float zeroTwo = zeroTerm * twoTerm / denominatorTerm;
                    float oneOne = oneTerm * oneTerm / denominatorTerm;
                    float oneTwo = oneTerm * twoTerm / denominatorTerm;
                    float twoTwo = twoTerm * twoTerm / denominatorTerm;

                    for (int a = 0; a < 4; a++) {
                        for (int b = 0; b < 4; b++) {
                            for (int c = 0; c < 4; c++) {
                                for (int d = 0; d < 4; d++) {
                                    int temp = a << 6 | b << 4 | c << 2 | d;
                                    answer[temp] = zeroZero;
                                    answer[0x100 | temp] = zeroOne;
                                    answer[0x400 | temp] = zeroOne;
                                    answer[0x200 | temp] = zeroTwo;
                                    answer[0x800 | temp] = zeroTwo;
                                    answer[0x500 | temp] = oneOne;
                                    answer[0x600 | temp] = oneTwo;
                                    answer[0x900 | temp] = oneTwo;
                                    answer[0xA00 | temp] = twoTwo;
                                }
                            }
                        }
                    }

                    for (int i = 0; i < myNumTaxa; i++) {
                        byte genotype = myGenotypes.genotype(i, myCurrentSite);
                        majorCount[i] = (short) (PRECALCULATED_COUNTS[((major & 0x7) << 6) | ((genotype & 0x70) >>> 1) | (genotype & 0x7)] << 8);
                    }
                } else {
                    for (int t = 0; t < myNumTaxa; t++) {
                        majorCount[t] = 0x300;
                    }
                }

                if (numSitesPerBlock > 1) {

                    int currentSite = myCurrentSite + 1;
                    major = myGenotypes.majorAllele(currentSite);
                    majorFreq = (float) myGenotypes.majorAlleleFrequency(currentSite);
                    majorFreqTimes2 = majorFreq * 2.0f;
                    denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                    if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                        float zeroTerm = 0.0f - majorFreqTimes2;
                        float oneTerm = 1.0f - majorFreqTimes2;
                        float twoTerm = 2.0f - majorFreqTimes2;

                        float zeroZero = zeroTerm * zeroTerm / denominatorTerm;
                        float zeroOne = zeroTerm * oneTerm / denominatorTerm;
                        float zeroTwo = zeroTerm * twoTerm / denominatorTerm;
                        float oneOne = oneTerm * oneTerm / denominatorTerm;
                        float oneTwo = oneTerm * twoTerm / denominatorTerm;
                        float twoTwo = twoTerm * twoTerm / denominatorTerm;

                        for (int a = 0; a < 4; a++) {
                            for (int b = 0; b < 4; b++) {
                                for (int c = 0; c < 4; c++) {
                                    for (int d = 0; d < 4; d++) {
                                        int temp = a << 10 | b << 8 | c << 2 | d;
                                        answer[temp] += zeroZero;
                                        answer[temp | 0x10] += zeroOne;
                                        answer[temp | 0x40] += zeroOne;
                                        answer[temp | 0x20] += zeroTwo;
                                        answer[temp | 0x80] += zeroTwo;
                                        answer[temp | 0x50] += oneOne;
                                        answer[temp | 0x60] += oneTwo;
                                        answer[temp | 0x90] += oneTwo;
                                        answer[temp | 0xA0] += twoTwo;
                                    }
                                }
                            }
                        }

                        for (int i = 0; i < myNumTaxa; i++) {
                            byte genotype = myGenotypes.genotype(i, currentSite);
                            majorCount[i] = (short) (majorCount[i] | PRECALCULATED_COUNTS[((major & 0x7) << 6) | ((genotype & 0x70) >>> 1) | (genotype & 0x7)] << 4);
                        }
                    } else {
                        for (int t = 0; t < myNumTaxa; t++) {
                            majorCount[t] |= 0x30;
                        }
                    }
                } else {
                    for (int t = 0; t < myNumTaxa; t++) {
                        majorCount[t] |= 0x30;
                    }
                }

                if (numSitesPerBlock > 2) {

                    int currentSite = myCurrentSite + 2;
                    major = myGenotypes.majorAllele(currentSite);
                    majorFreq = (float) myGenotypes.majorAlleleFrequency(currentSite);
                    majorFreqTimes2 = majorFreq * 2.0f;
                    denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                    if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                        float zeroTerm = 0.0f - majorFreqTimes2;
                        float oneTerm = 1.0f - majorFreqTimes2;
                        float twoTerm = 2.0f - majorFreqTimes2;

                        float zeroZero = zeroTerm * zeroTerm / denominatorTerm;
                        float zeroOne = zeroTerm * oneTerm / denominatorTerm;
                        float zeroTwo = zeroTerm * twoTerm / denominatorTerm;
                        float oneOne = oneTerm * oneTerm / denominatorTerm;
                        float oneTwo = oneTerm * twoTerm / denominatorTerm;
                        float twoTwo = twoTerm * twoTerm / denominatorTerm;

                        for (int a = 0; a < 4; a++) {
                            for (int b = 0; b < 4; b++) {
                                for (int c = 0; c < 4; c++) {
                                    for (int d = 0; d < 4; d++) {
                                        int temp = a << 10 | b << 8 | c << 6 | d << 4;
                                        answer[temp] += zeroZero;
                                        answer[temp | 0x1] += zeroOne;
                                        answer[temp | 0x4] += zeroOne;
                                        answer[temp | 0x2] += zeroTwo;
                                        answer[temp | 0x8] += zeroTwo;
                                        answer[temp | 0x5] += oneOne;
                                        answer[temp | 0x6] += oneTwo;
                                        answer[temp | 0x9] += oneTwo;
                                        answer[temp | 0xA] += twoTwo;
                                    }
                                }
                            }
                        }

                        for (int i = 0; i < myNumTaxa; i++) {
                            byte genotype = myGenotypes.genotype(i, currentSite);
                            majorCount[i] = (short) (majorCount[i] | PRECALCULATED_COUNTS[((major & 0x7) << 6) | ((genotype & 0x70) >>> 1) | (genotype & 0x7)]);
                        }
                    } else {
                        for (int t = 0; t < myNumTaxa; t++) {
                            majorCount[t] |= 0x3;
                        }
                    }
                } else {
                    for (int t = 0; t < myNumTaxa; t++) {
                        majorCount[t] |= 0x3;
                    }
                }

                int index = 0;
                for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                    if (majorCount[firstTaxa] != 0x333) {
                        int temp = majorCount[firstTaxa] << 2;
                        for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                            int aIndex = temp | majorCount[secondTaxa];
                            distances[index] += answer[aIndex];
                            counts[index] += INCREMENT[aIndex];
                            index++;
                        }
                    } else {
                        index += myNumTaxa - firstTaxa;
                    }
                }
            }

            action.accept(result);
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);
        }

        @Override
        public boolean tryAdvance(Consumer<? super CountersDistances> action) {
            if (myCurrentSite < myFence) {

                CountersDistances result = new CountersDistances(myNumTaxa);
                int[] counts = result.myCounters;
                float[] distances = result.myDistances;
                byte[] majorCount = new byte[myNumTaxa];
                float[] answer = new float[12];
                byte[] increment = new byte[12];
                increment[0] = increment[1] = increment[2]
                        = increment[4] = increment[5] = increment[6]
                        = increment[8] = increment[9] = increment[10] = 1;

                byte major = myGenotypes.majorAllele(myCurrentSite);
                float majorFreq = (float) myGenotypes.majorAlleleFrequency(myCurrentSite);
                float majorFreqTimes2 = majorFreq * 2.0f;
                float denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);
                if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                    float zeroTerm = 0.0f - majorFreqTimes2;
                    float oneTerm = 1.0f - majorFreqTimes2;
                    float twoTerm = 2.0f - majorFreqTimes2;

                    answer[0] = zeroTerm * zeroTerm / denominatorTerm;
                    answer[1] = answer[4] = zeroTerm * oneTerm / denominatorTerm;
                    answer[2] = answer[8] = zeroTerm * twoTerm / denominatorTerm;
                    answer[5] = oneTerm * oneTerm / denominatorTerm;
                    answer[6] = answer[9] = oneTerm * twoTerm / denominatorTerm;
                    answer[10] = twoTerm * twoTerm / denominatorTerm;

                    for (int i = 0; i < myNumTaxa; i++) {
                        byte genotype = myGenotypes.genotype(i, myCurrentSite);
                        if (genotype == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                            majorCount[i] = 3;
                        } else {
                            majorCount[i] = 0;
                            if ((genotype & 0xF) == major) {
                                majorCount[i]++;
                            }
                            if (((genotype >>> 4) & 0xF) == major) {
                                majorCount[i]++;
                            }
                        }
                    }

                    int index = 0;
                    for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                        if (majorCount[firstTaxa] != 3) {
                            int temp = majorCount[firstTaxa] << 2;
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                int aIndex = temp | majorCount[secondTaxa];
                                distances[index] += answer[aIndex];
                                counts[index] += increment[aIndex];
                                index++;
                            }
                        } else {
                            index += myNumTaxa - firstTaxa;
                        }
                    }
                }

                action.accept(result);

                return true;
            } else {
                return false;
            }
        }

        @Override
        public Spliterator<CountersDistances> trySplit() {
            int lo = myCurrentSite;
            int mid = (lo + myFence) >>> 1;
            if (lo < mid) {
                myCurrentSite = mid;
                return new GCTASiteSpliterator(myGenotypes, lo, mid, myProgressListener);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return (long) (myFence - myCurrentSite);
        }

        @Override
        public int characteristics() {
            return IMMUTABLE;
        }
    }

}
