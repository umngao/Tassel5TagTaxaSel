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
        int[][] counts = counters.myCounters;
        double[][] distances = counters.myDistances;

        double[][] result = new double[numSeqs][numSeqs];

        for (int t = 0; t < numSeqs; t++) {
            for (int i = 0, n = counts[t].length; i < n; i++) {
                result[t][t + i] = result[t + i][t] = distances[t][i] / (double) counts[t][i];
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

        private final int[][] myCounters;
        private final double[][] myDistances;
        private final int myNumTaxa;

        public CountersDistances(int numTaxa) {
            myNumTaxa = numTaxa;
            myCounters = new int[myNumTaxa][];
            myDistances = new double[myNumTaxa][];
            for (int i = 0; i < myNumTaxa; i++) {
                myCounters[i] = new int[(myNumTaxa - i)];
                myDistances[i] = new double[(myNumTaxa - i)];
            }
        }

        public void addAll(CountersDistances counters) {
            double[][] otherDistances = counters.myDistances;
            for (int t = 0; t < myNumTaxa; t++) {
                for (int i = 0, n = myNumTaxa - t; i < n; i++) {
                    myDistances[t][i] += otherDistances[t][i];
                }
            }
            otherDistances = null;
            int[][] otherCounters = counters.myCounters;
            for (int t = 0; t < myNumTaxa; t++) {
                for (int i = 0, n = myNumTaxa - t; i < n; i++) {
                    myCounters[t][i] += otherCounters[t][i];
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
            int[][] counts = result.myCounters;
            double[][] distances = result.myDistances;
            double[] term = new double[myNumTaxa];
            byte[] countTerm = new byte[myNumTaxa];
            for (; myCurrentSite < myFence; myCurrentSite++) {

                byte major = myGenotypes.majorAllele(myCurrentSite);
                double majorFreq = myGenotypes.majorAlleleFrequency(myCurrentSite);
                double majorFreqTimes2 = majorFreq * 2.0;
                double denominatorTerm = majorFreqTimes2 * (1.0 - majorFreq);
                if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                    byte majorTimes2 = (byte) ((major << 4) | major);
                    double zeroTerm = 0.0 - majorFreqTimes2;
                    double oneTerm = 1.0 - majorFreqTimes2;
                    double twoTerm = 2.0 - majorFreqTimes2;

                    for (int i = 0; i < myNumTaxa; i++) {
                        byte genotype = myGenotypes.genotype(i, myCurrentSite);
                        if (genotype == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                            term[i] = 0.0;
                            countTerm[i] = 0;
                        } else if (genotype == majorTimes2) {
                            term[i] = twoTerm;
                            countTerm[i] = 1;
                        } else if ((genotype & 0xF) == major) {
                            term[i] = oneTerm;
                            countTerm[i] = 1;
                        } else if (((genotype >>> 4) & 0xF) == major) {
                            term[i] = oneTerm;
                            countTerm[i] = 1;
                        } else {
                            term[i] = zeroTerm;
                            countTerm[i] = 1;
                        }
                    }

                    for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                        if (countTerm[firstTaxa] != 0) {
                            double firstTerm = term[firstTaxa] / denominatorTerm;
                            int index = 0;
                            double[] distanceRow = distances[firstTaxa];
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                distanceRow[index] += (firstTerm * term[secondTaxa]);
                                index++;
                            }
                            index = 0;
                            int[] countRow = counts[firstTaxa];
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                countRow[index] += countTerm[secondTaxa];
                                index++;
                            }
                        }
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
                int[][] counts = result.myCounters;
                double[][] distances = result.myDistances;
                double[] term = new double[myNumTaxa];
                byte[] countTerm = new byte[myNumTaxa];

                byte major = myGenotypes.majorAllele(myCurrentSite);
                double majorFreq = myGenotypes.majorAlleleFrequency(myCurrentSite);
                double majorFreqTimes2 = majorFreq * 2.0;
                double denominatorTerm = majorFreqTimes2 * (1.0 - majorFreq);
                if ((major != GenotypeTable.UNKNOWN_ALLELE) && (denominatorTerm != 0.0)) {

                    byte majorTimes2 = (byte) ((major << 4) | major);
                    double zeroTerm = 0.0 - majorFreqTimes2;
                    double oneTerm = 1.0 - majorFreqTimes2;
                    double twoTerm = 2.0 - majorFreqTimes2;

                    for (int i = 0; i < myNumTaxa; i++) {
                        byte genotype = myGenotypes.genotype(i, myCurrentSite);
                        if (genotype == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                            term[i] = 0.0;
                            countTerm[i] = 0;
                        } else if (genotype == majorTimes2) {
                            term[i] = twoTerm;
                            countTerm[i] = 1;
                        } else if ((genotype & 0xF) == major) {
                            term[i] = oneTerm;
                            countTerm[i] = 1;
                        } else if (((genotype >>> 4) & 0xF) == major) {
                            term[i] = oneTerm;
                            countTerm[i] = 1;
                        } else {
                            term[i] = zeroTerm;
                            countTerm[i] = 1;
                        }
                    }

                    for (int firstTaxa = 0; firstTaxa < myNumTaxa; firstTaxa++) {
                        if (countTerm[firstTaxa] != 0) {
                            double firstTerm = term[firstTaxa] / denominatorTerm;
                            int index = 0;
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                distances[firstTaxa][index] += (firstTerm * term[secondTaxa]);
                                index++;
                            }
                            index = 0;
                            for (int secondTaxa = firstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                counts[firstTaxa][index] += countTerm[secondTaxa];
                                index++;
                            }
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
