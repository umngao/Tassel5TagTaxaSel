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
            int[][] otherCounters = counters.myCounters;
            double[][] otherDistances = counters.myDistances;
            for (int t = 0; t < myNumTaxa; t++) {
                for (int i = 0, n = myCounters[t].length; i < n; i++) {
                    myCounters[t][i] += otherCounters[t][i];
                    myDistances[t][i] += otherDistances[t][i];
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
        private int myFirstTaxa;
        private final ProgressListener myProgressListener;

        GCTASiteSpliterator(GenotypeTable genotypes, int currentIndex, int fence, ProgressListener listener) {
            myGenotypes = genotypes;
            myNumTaxa = myGenotypes.numberOfTaxa();
            myNumSites = myGenotypes.numberOfSites();
            myFirstTaxa = 0;
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

                    for (int i = 0; i < myNumTaxa; i++) {
                        byte genotype = myGenotypes.genotype(i, myCurrentSite);
                        int key = (genotype & 0x70) >>> 1 | (genotype & 0x7);
                        if (key == 0x3F) {
                            term[i] = 0.0;
                            countTerm[i] = 0;
                        } else {
                            byte count = 0;
                            if ((genotype & 0xF) == major) {
                                count++;
                            }
                            if ((genotype >>> 4) == major) {
                                count++;
                            }
                            term[i] = (double) count - majorFreqTimes2;
                            countTerm[i] = 1;
                        }
                    }

                    for (; myFirstTaxa < myNumTaxa;) {
                        if (countTerm[myFirstTaxa] != 0) {
                            double firstTerm = term[myFirstTaxa] / denominatorTerm;
                            int index = 0;
                            for (int secondTaxa = myFirstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                                counts[myFirstTaxa][index] += countTerm[secondTaxa];
                                distances[myFirstTaxa][index] += (firstTerm * term[secondTaxa]);
                                index++;
                            }
                        }
                        myFirstTaxa++;
                    }
                    myFirstTaxa = 0;
                }
            }
            action.accept(result);
            myNumSitesProcessed += numSitesProcessed;
            fireProgress((int) ((double) myNumSitesProcessed / (double) myNumSites * 100.0), myProgressListener);
        }

        @Override
        public boolean tryAdvance(Consumer<? super CountersDistances> action) {
            if ((myCurrentSite < myFence) && (myFirstTaxa < myNumTaxa)) {
                //int key1 = (myCachedSiteGenotype[myFirstTaxa] & 0x70) >>> 1 | (myCachedSiteGenotype[myFirstTaxa] & 0x7);
                //long[] incrementFunctions = INCREMENT_FUNCTIONS[key1];

                long[] result = new long[myNumTaxa - myFirstTaxa];
                int index = 0;
                for (int secondTaxa = myFirstTaxa; secondTaxa < myNumTaxa; secondTaxa++) {
                    //int key2 = (myCachedSiteGenotype[secondTaxa] & 0x70) >>> 1 | (myCachedSiteGenotype[secondTaxa] & 0x7);
                    //result[index++] = incrementFunctions[key2];
                }
                //action.accept(result);

                myFirstTaxa++;
                if (myFirstTaxa == myNumTaxa) {
                    myCurrentSite++;
                    if (myCurrentSite < myFence) {
                        //myCachedSiteGenotype = myGenotypes.genotypeAllTaxa(myCurrentSite);
                        //myCachedSite = myCurrentSite;
                    } else {
                        //myCachedSiteGenotype = null;
                        //myCachedSite = -1;
                    }
                    myFirstTaxa = 0;
                }
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
