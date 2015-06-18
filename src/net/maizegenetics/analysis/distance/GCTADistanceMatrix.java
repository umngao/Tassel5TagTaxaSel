/*
 *  GCTADistanceMatrix
 * 
 *  Created on May 31, 2015
 */
package net.maizegenetics.analysis.distance;

import java.util.Arrays;
import java.util.Optional;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Tuple;

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

    private static final byte[] INCREMENT = new byte[65536];

    static {
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                int temp = a << 14 | b << 12;
                for (int c = 0; c < 4; c++) {
                    for (int d = 0; d < 4; d++) {
                        int temp2 = c << 10 | d << 8;
                        for (int e = 0; e < 4; e++) {
                            for (int f = 0; f < 4; f++) {
                                int temp3 = e << 6 | f << 4;
                                for (int g = 0; g < 4; g++) {
                                    for (int h = 0; h < 4; h++) {
                                        int incrementIndex = temp | temp2 | temp3 | g << 2 | h;
                                        if ((a != 3) && (b != 3)) {
                                            INCREMENT[incrementIndex]++;
                                        }
                                        if ((c != 3) && (d != 3)) {
                                            INCREMENT[incrementIndex]++;
                                        }
                                        if ((e != 3) && (f != 3)) {
                                            INCREMENT[incrementIndex]++;
                                        }
                                        if ((g != 3) && (h != 3)) {
                                            INCREMENT[incrementIndex]++;
                                        }
                                    }
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

            float[] answer1 = new float[65536];
            float[] answer2 = new float[65536];
            float[] answer3 = new float[65536];

            for (; myCurrentSite < myFence; myCurrentSite += 12) {

                Tuple<short[], float[]> firstThree = getThreeSites(myCurrentSite);
                float[] possibleTerms = firstThree.y;
                short[] majorCount1 = firstThree.x;

                Tuple<short[], float[]> secondThree = getThreeSites(myCurrentSite + 4);
                float[] possibleTerms2 = secondThree.y;
                short[] majorCount2 = secondThree.x;

                Tuple<short[], float[]> thirdThree = getThreeSites(myCurrentSite + 8);
                float[] possibleTerms3 = thirdThree.y;
                short[] majorCount3 = thirdThree.x;

                for (int i = 0; i < 65536; i++) {
                    answer1[i] = possibleTerms[(i & 0xF000) >>> 12] + possibleTerms[((i & 0xF00) >>> 8) | 0x10] + possibleTerms[((i & 0xF0) >>> 4) | 0x20] + possibleTerms[(i & 0xF) | 0x30];
                    answer2[i] = possibleTerms2[(i & 0xF000) >>> 12] + possibleTerms2[((i & 0xF00) >>> 8) | 0x10] + possibleTerms2[((i & 0xF0) >>> 4) | 0x20] + possibleTerms2[(i & 0xF) | 0x30];
                    answer3[i] = possibleTerms3[(i & 0xF000) >>> 12] + possibleTerms3[((i & 0xF00) >>> 8) | 0x10] + possibleTerms3[((i & 0xF0) >>> 4) | 0x20] + possibleTerms3[(i & 0xF) | 0x30];
                }

                int index = 0;
                for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                    if ((majorCount1[firstTaxa] != 0x3333) || (majorCount2[firstTaxa] != 0x3333) || (majorCount3[firstTaxa] != 0x3333)) {
                        int temp1 = majorCount1[firstTaxa] << 2;
                        int temp2 = majorCount2[firstTaxa] << 2;
                        int temp3 = majorCount3[firstTaxa] << 2;
                        for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                            int aIndex = temp1 | majorCount1[secondTaxa];
                            int bIndex = temp2 | majorCount2[secondTaxa];
                            int cIndex = temp3 | majorCount3[secondTaxa];
                            distances[index] += answer1[aIndex] + answer2[bIndex] + answer3[cIndex];
                            counts[index] += INCREMENT[aIndex] + INCREMENT[bIndex] + INCREMENT[cIndex];
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

        private Tuple<short[], float[]> getThreeSites(int currentSite) {

            float[] possibleTerms = new float[64];
            short[] majorCount = new short[myNumTaxa];
            Arrays.fill(majorCount, (short) 0x3333);

            if (currentSite < myFence) {

                byte major = myGenotypes.majorAllele(currentSite);
                float majorFreq = (float) myGenotypes.majorAlleleFrequency(currentSite);
                float majorFreqTimes2 = majorFreq * 2.0f;
                float denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                float[] term = new float[3];

                if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                    term[0] = 0.0f - majorFreqTimes2;
                    term[1] = 1.0f - majorFreqTimes2;
                    term[2] = 2.0f - majorFreqTimes2;

                    possibleTerms[0] = term[0] * term[0] / denominatorTerm;
                    possibleTerms[1] = possibleTerms[4] = term[0] * term[1] / denominatorTerm;
                    possibleTerms[2] = possibleTerms[8] = term[0] * term[2] / denominatorTerm;
                    possibleTerms[5] = term[1] * term[1] / denominatorTerm;
                    possibleTerms[6] = possibleTerms[9] = term[1] * term[2] / denominatorTerm;
                    possibleTerms[10] = term[2] * term[2] / denominatorTerm;

                    int temp = (major & 0x7) << 6;
                    for (int i = 0; i < myNumTaxa; i++) {
                        byte genotype = myGenotypes.genotype(i, currentSite);
                        majorCount[i] = (short) (0x333 | PRECALCULATED_COUNTS[temp | ((genotype & 0x70) >>> 1) | (genotype & 0x7)] << 12);
                    }
                }

                currentSite++;
                if (currentSite < myFence) {

                    major = myGenotypes.majorAllele(currentSite);
                    majorFreq = (float) myGenotypes.majorAlleleFrequency(currentSite);
                    majorFreqTimes2 = majorFreq * 2.0f;
                    denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                    if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                        term[0] = 0.0f - majorFreqTimes2;
                        term[1] = 1.0f - majorFreqTimes2;
                        term[2] = 2.0f - majorFreqTimes2;

                        possibleTerms[16] = term[0] * term[0] / denominatorTerm;
                        possibleTerms[17] = possibleTerms[20] = term[0] * term[1] / denominatorTerm;
                        possibleTerms[18] = possibleTerms[24] = term[0] * term[2] / denominatorTerm;
                        possibleTerms[21] = term[1] * term[1] / denominatorTerm;
                        possibleTerms[22] = possibleTerms[25] = term[1] * term[2] / denominatorTerm;
                        possibleTerms[26] = term[2] * term[2] / denominatorTerm;

                        int temp = (major & 0x7) << 6;
                        for (int i = 0; i < myNumTaxa; i++) {
                            byte genotype = myGenotypes.genotype(i, currentSite);
                            majorCount[i] = (short) (majorCount[i] & (0x3033 | PRECALCULATED_COUNTS[temp | ((genotype & 0x70) >>> 1) | (genotype & 0x7)] << 8));
                        }
                    }

                    currentSite++;
                    if (currentSite < myFence) {

                        major = myGenotypes.majorAllele(currentSite);
                        majorFreq = (float) myGenotypes.majorAlleleFrequency(currentSite);
                        majorFreqTimes2 = majorFreq * 2.0f;
                        denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                        if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                            term[0] = 0.0f - majorFreqTimes2;
                            term[1] = 1.0f - majorFreqTimes2;
                            term[2] = 2.0f - majorFreqTimes2;

                            possibleTerms[32] = term[0] * term[0] / denominatorTerm;
                            possibleTerms[33] = possibleTerms[36] = term[0] * term[1] / denominatorTerm;
                            possibleTerms[34] = possibleTerms[40] = term[0] * term[2] / denominatorTerm;
                            possibleTerms[37] = term[1] * term[1] / denominatorTerm;
                            possibleTerms[38] = possibleTerms[41] = term[1] * term[2] / denominatorTerm;
                            possibleTerms[42] = term[2] * term[2] / denominatorTerm;

                            int temp = (major & 0x7) << 6;
                            for (int i = 0; i < myNumTaxa; i++) {
                                byte genotype = myGenotypes.genotype(i, currentSite);
                                majorCount[i] = (short) (majorCount[i] & (0x3303 | PRECALCULATED_COUNTS[temp | ((genotype & 0x70) >>> 1) | (genotype & 0x7)] << 4));
                            }
                        }

                        currentSite++;
                        if (currentSite < myFence) {

                            major = myGenotypes.majorAllele(currentSite);
                            majorFreq = (float) myGenotypes.majorAlleleFrequency(currentSite);
                            majorFreqTimes2 = majorFreq * 2.0f;
                            denominatorTerm = majorFreqTimes2 * (1.0f - majorFreq);

                            if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                                term[0] = 0.0f - majorFreqTimes2;
                                term[1] = 1.0f - majorFreqTimes2;
                                term[2] = 2.0f - majorFreqTimes2;

                                possibleTerms[48] = term[0] * term[0] / denominatorTerm;
                                possibleTerms[49] = possibleTerms[52] = term[0] * term[1] / denominatorTerm;
                                possibleTerms[50] = possibleTerms[56] = term[0] * term[2] / denominatorTerm;
                                possibleTerms[53] = term[1] * term[1] / denominatorTerm;
                                possibleTerms[54] = possibleTerms[57] = term[1] * term[2] / denominatorTerm;
                                possibleTerms[58] = term[2] * term[2] / denominatorTerm;

                                int temp = (major & 0x7) << 6;
                                for (int i = 0; i < myNumTaxa; i++) {
                                    byte genotype = myGenotypes.genotype(i, currentSite);
                                    majorCount[i] = (short) (majorCount[i] & (0x3330 | PRECALCULATED_COUNTS[temp | ((genotype & 0x70) >>> 1) | (genotype & 0x7)]));
                                }
                            }
                        }
                    }

                }

            }

            return new Tuple<>(majorCount, possibleTerms);

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
