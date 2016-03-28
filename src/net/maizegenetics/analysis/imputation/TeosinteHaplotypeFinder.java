package net.maizegenetics.analysis.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;

import net.maizegenetics.analysis.clustering.Haplotype;
import net.maizegenetics.analysis.clustering.HaplotypeClusterer;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;

public class TeosinteHaplotypeFinder extends BiparentalHaplotypeFinder {
    private PopulationData myPopulationData;

    public TeosinteHaplotypeFinder(PopulationData popdata) {
        super(popdata);
        myPopulationData = popdata;
    }
    
    public void assignHaplotyes() {
        //test diagnostic output
//        List<Position> plist = new ArrayList<>();
//        for (int i = 0; i < window; i++) {
//                plist.add(new GeneralPosition.Builder(new Chromosome("0"), i).build());
//        }
        
        int startIncr = window - overlap;
        int diff = maxDifferenceScore;
        
        //assign haplotype of parent1(A), parent2(C), het(M) at each non-missing locus for each taxon
        //for each window:
        //      1. cluster sequence
        //      2. sort clusters by score
        //      3. move haplotypes to largest consistent cluster
        //      4. assign all clusters > min size to one of parents
        //  5. assign remaining clusters to parent or het
        //      6. record parent or het for each taxon and non-missing site in the window
        
        GenotypeTable filterGeno = preFilterSites();
        
        myPopulationData.original = filterGeno;
        
        int nsites = myPopulationData.original.numberOfSites();
        myPopulationData.alleleA = new byte[nsites];
        myPopulationData.alleleC = new byte[nsites];
        Arrays.fill(myPopulationData.alleleA, NN);
        Arrays.fill(myPopulationData.alleleC, NN);
        
        //find a window with exactly two haplotypes
        boolean exactlyTwo = false;
        int initialStart = 0;
        int maxStart = nsites - window;
        
        //this structure allows a parent to have more than one haplotype in a segment
        ArrayList<ArrayList<Haplotype>> parentHaplotypes = new ArrayList<ArrayList<Haplotype>>(2);  
        parentHaplotypes.add(new ArrayList<Haplotype>());
        parentHaplotypes.add(new ArrayList<Haplotype>());
        Haplotype h0 = null;
        Haplotype h1 = null;
        while (!exactlyTwo & initialStart < maxStart) {
                int minNotMissing = (int) (window * minNotMissingProportion);
                HaplotypeClusterer myClusterMaker = clusterWindow(filterGeno, initialStart, window, diff, minNotMissing); //1
                myClusterMaker.sortClusters(); //2
                myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
                myClusterMaker.removeHeterozygousClusters(5 + diff);
                
                h0 = new Haplotype(myClusterMaker.getClusterList().get(0).getMajorityHaplotype());
                h1 = new Haplotype(myClusterMaker.getClusterList().get(1).getMajorityHaplotype());
                int cluster2Size = 0;
                if (myClusterMaker.getClusterList().size() > 2) cluster2Size = myClusterMaker.getClusterList().get(2).getSize();
                if (cluster2Size < minClusterSize && h0.distanceFrom(h1) >= 2 * window - 4) {
                        exactlyTwo = true;
                        parentHaplotypes.get(0).add(h0);
                        parentHaplotypes.get(1).add(h1);
                        addHaplotypesToReport(myClusterMaker, initialStart);
                        
                } else {
                        initialStart += window;
                }
        }
        
        if (!exactlyTwo) {
                throw new RuntimeException("Unable to find start window with only two haplotypes.");
        }
        
        //update parent haplotype alleles in myPopulationData
        ArrayList<Position> filterPositions = new ArrayList<Position>();
        for (int s = initialStart; s < initialStart + window; s++) filterPositions.add(filterGeno.positions().get(s));
        updatePopulationDataAlleles(parentHaplotypes, filterPositions, 0, window);
        
        for (int start = initialStart + startIncr; start < nsites - overlap; start += startIncr) {
                int windowSize = window;
                if (start + window >= nsites) {
                    windowSize = nsites - start;
                    diff = diff * windowSize / window;
                }
                int minNotMissingAdjusted = (int) (windowSize * minNotMissingProportion);
                HaplotypeClusterer myClusterMaker = clusterWindow(filterGeno, start, windowSize, diff, minNotMissingAdjusted); //1

                myClusterMaker.sortClusters(); //2
                myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
                myClusterMaker.removeHeterozygousClusters(diff + 5);
                addHaplotypesToReport(myClusterMaker, start);

                //get major haplotypes
                ArrayList<Haplotype> myHaplotypes = mergeMajorHaplotypes(myClusterMaker, minClusterSize);
                
                //assign parent to each Haplotype based on previous Haplotypes
                parentHaplotypes = getParentHaplotypes(parentHaplotypes, myHaplotypes, overlap, true);
                
                //update parent haplotype alleles in myPopulationData
                filterPositions.clear();
                for (int s = start; s < start + windowSize; s++) filterPositions.add(filterGeno.positions().get(s));
                updatePopulationDataAlleles(parentHaplotypes, filterPositions, overlap, windowSize - overlap);
        }
        
        //do the same thing starting at the original window going in the reverse direction
        parentHaplotypes.get(0).clear();
        parentHaplotypes.get(0).add(h0);
        parentHaplotypes.get(1).clear();
        parentHaplotypes.get(1).add(h1);
        diff = maxDifferenceScore;
        for (int start = (initialStart - startIncr); start > -startIncr; start -= startIncr) {
                int end = start + window;
                if (start < 0) {
                    start = 0;
                    diff = diff * (end - start) / window;
                }
                int windowSize = end - start;
                int minNotMissingAdjusted = (int) (windowSize * minNotMissingProportion);
                HaplotypeClusterer myClusterMaker = clusterWindow(myPopulationData.original, start, windowSize, diff, minNotMissingAdjusted); //1
                myClusterMaker.sortClusters(); //2
                myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
                myClusterMaker.removeHeterozygousClusters(5 + diff);
                addHaplotypesToReport(myClusterMaker, start);
                
                //get major haplotypes
                ArrayList<Haplotype> myHaplotypes = mergeMajorHaplotypes(myClusterMaker, minClusterSize);
                
                //assign parent to each Haplotype based on previous Haplotypes
                parentHaplotypes = getParentHaplotypes(parentHaplotypes, myHaplotypes, overlap, false);

                //update parent haplotype alleles in myPopulationData
                filterPositions.clear();
                for (int s = start; s < start + windowSize; s++) filterPositions.add(filterGeno.positions().get(s));
                updatePopulationDataAlleles(parentHaplotypes, filterPositions, 0, windowSize - overlap);
        }
        
}

    @Override
    public GenotypeTable preFilterSites() {
        int ntaxa = myPopulationData.original.numberOfTaxa();
        GenotypeTable filterGeno = NucleotideImputationUtils.filterSnpsByTag(myPopulationData.original, minMaf, 1 - minCoverage, 1);
        
        
        return filterGeno;
    }

    
}
