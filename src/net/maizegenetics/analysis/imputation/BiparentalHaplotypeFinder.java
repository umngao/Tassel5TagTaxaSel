package net.maizegenetics.analysis.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeSet;

import org.apache.commons.math.stat.StatUtils;

import net.maizegenetics.analysis.clustering.Haplotype;
import net.maizegenetics.analysis.clustering.HaplotypeCluster;
import net.maizegenetics.analysis.clustering.HaplotypeClusterer;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 * @author pbradbury
 * 
 */
public class BiparentalHaplotypeFinder {
	PopulationData myPopulationData;
	static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	static final byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
	static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	
	//parameters
	/**
	 * The size of the window in which to evaluate haplotypes
	 */
	int window = 100;
	
	/**
	 * The extent of overlap between adjacent windows. The overlap is used to assign individual haplotypes to the same parent groups as the previous window.
	 */
	int overlap = 25;
	
	/**
	 * Only cluster haplotypes with a minimum number of non missing values.
	 */
	int minNotMissing = 20;
	
	/**
	 * Only use clusters with minClusterSize taxa as parent haplotypes
	 */
	int minClusterSize = 5;
	
	/**
	 * Filter out sites with less than minR2 average r-square with neighboring sites (window size = 50).
	 */
	double minR2 = 0.2;
	
	/**
	 * Filter out sites with minimum allele frequency less than minMaf.
	 */
	double minMaf = 0.05;
	
	/**
	 * Filter out sites with coverage less than minCoverage.
	 */
	double minCoverage = 0.2;
	
	/**
	 * Filter out sites more than maxHetDeviation * (standard deviation) different from the mean percent heterozygosity.;
	 */
	double maxHetDeviation = 5;
	
	public BiparentalHaplotypeFinder(PopulationData popdata) {
		myPopulationData = popdata;
		
	}
	
	public void assignHaplotyes() {
		int startIncr = window - overlap;
		int diff = 0;
		
		//assign haplotype of parent1(A), parent2(C), het(M) at each non-missing locus for each taxon
		//for each window:
		//	1. cluster sequence
		//	2. sort clusters by score
		//	3. move haplotypes to largest consistent cluster
		//	4. assign all clusters > min size to one of parents
		//  5. assign remaining clusters to parent or het
		//	6. record parent or het for each taxon and non-missing site in the window
		
//		GenotypeTable filterGeno = myPopulationData.original; 
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
			HaplotypeClusterer myClusterMaker = clusterWindow(filterGeno, initialStart, window, diff, minNotMissing); //1
			myClusterMaker.sortClusters(); //2
			myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
			h0 = new Haplotype(myClusterMaker.getClusterList().get(0).getMajorityHaplotype());
			h1 = new Haplotype(myClusterMaker.getClusterList().get(1).getMajorityHaplotype());
			int cluster2Size = myClusterMaker.getClusterList().get(2).getSize();
			if (cluster2Size < minClusterSize && h0.distanceFrom(h1) >= 2 * window - 4) {
				exactlyTwo = true;
				parentHaplotypes.get(0).add(h0);
				parentHaplotypes.get(1).add(h1);
				
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
			if (start + window >= nsites) windowSize = nsites - start;
			int minNotMissingAdjusted = (int) Math.floor(minNotMissing * ((double) windowSize)/window);
			HaplotypeClusterer myClusterMaker = clusterWindow(filterGeno, start, windowSize, diff, minNotMissingAdjusted); //1
			
			//debug
			myClusterMaker.sortClusters(); //2
			myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
			
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
		for (int start = (initialStart - startIncr); start > -startIncr; start -= startIncr) {
			int end = start + window;
			if (start < 0) start = 0;
			int windowSize = end - start;
			HaplotypeClusterer myClusterMaker = clusterWindow(myPopulationData.original, start, windowSize, diff, minNotMissing); //1
			myClusterMaker.sortClusters(); //2
			myClusterMaker.moveAllHaplotypesToBiggestCluster(diff); //3
			
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
	
	public HaplotypeClusterer clusterWindow(GenotypeTable a, int start, int length, int maxdif, int minNotMissingPerHaplotype) {
		int ntaxa = a.numberOfTaxa();
		int end = start + length;
		ArrayList<Haplotype> haps = new ArrayList<Haplotype>();
		
		for (int t = 0; t < ntaxa; t++) {
			Haplotype hap = new Haplotype(a.genotypeRange(t, start, end), t);
			if (hap.notMissingCount >= minNotMissingPerHaplotype) haps.add(hap);
		}
		
		
		HaplotypeClusterer hc = new HaplotypeClusterer(haps);
		hc.makeClusters();
		if (maxdif > 0) hc.mergeClusters(maxdif);
		return hc;
	}
	
	public void updatePopulationDataAlleles(ArrayList<ArrayList<Haplotype>> parentHaplotypes, ArrayList<Position> positions, int start, int length) {
	
		//at a site if one or both of the parents have multiple values, then any allele call not unique to a parent is set to N.
		int nhap0 = parentHaplotypes.get(0).size();
		int nhap1 = parentHaplotypes.get(1).size();
		if (nhap0 == 0 || nhap1 == 0) return;
		if (nhap0 == 1 && nhap1 == 1) {
			byte[] seqA = parentHaplotypes.get(0).get(0).seq;
			byte[] seqC = parentHaplotypes.get(1).get(0).seq;
			for (int ptr = start; ptr < start + length; ptr++) {
				int allelePtr = myPopulationData.original.positions().indexOf(positions.get(ptr));
				
				if (seqA[ptr] == seqC[ptr]) {
					myPopulationData.alleleA[allelePtr] = NN;
					myPopulationData.alleleC[allelePtr] = NN;
				} else {
					if (GenotypeTableUtils.isHeterozygous(seqA[ptr])) myPopulationData.alleleA[allelePtr] = NN;
					else myPopulationData.alleleA[allelePtr] = seqA[ptr];
					if (GenotypeTableUtils.isHeterozygous(seqC[ptr])) myPopulationData.alleleC[allelePtr] = NN;
					else myPopulationData.alleleC[allelePtr] = seqC[ptr];
				}
			}
		} else {
			for (int ptr = start; ptr < start + length; ptr++) {
				int allelePtr = myPopulationData.original.positions().indexOf(positions.get(ptr));
				TreeSet<Byte> p0alleles = new TreeSet<Byte>();
				TreeSet<Byte> p1alleles = new TreeSet<Byte>();
				for (Haplotype hap : parentHaplotypes.get(0)) {
					byte hapval = hap.seq[ptr];
					if (hapval != NN) {
						byte[] alleles = GenotypeTableUtils.getDiploidValues(hapval);
						p0alleles.add(alleles[0]);
						p0alleles.add(alleles[1]);
					}
				}
				for (Haplotype hap : parentHaplotypes.get(1)) {
					byte hapval = hap.seq[ptr];
					if (hapval != NN) {
						byte[] alleles = GenotypeTableUtils.getDiploidValues(hapval);
						p1alleles.add(alleles[0]);
						p1alleles.add(alleles[1]);
					}
				}
				
				if (p0alleles.size() == 0) { 
					myPopulationData.alleleA[allelePtr] = NN;
					if (p1alleles.size() == 0) myPopulationData.alleleC[allelePtr] = NN;
					else if (p1alleles.size() == 1) {
						byte allele = p1alleles.first();
						myPopulationData.alleleC[allelePtr] = GenotypeTableUtils.getDiploidValue(allele, allele);
					} else myPopulationData.alleleC[allelePtr] = NN;
				} else if (p0alleles.size() == 1) {
					byte Aallele = p0alleles.first();
					if (p1alleles.size() == 0) {
						myPopulationData.alleleA[allelePtr] = GenotypeTableUtils.getDiploidValue(Aallele, Aallele);
						myPopulationData.alleleC[allelePtr] = NN;
					} else if (p1alleles.size() == 1) {
						byte Callele = p1alleles.first();
						if (Aallele == Callele) {
							myPopulationData.alleleA[allelePtr] = NN;
							myPopulationData.alleleC[allelePtr] = NN;
						} else {
							myPopulationData.alleleA[allelePtr] = GenotypeTableUtils.getDiploidValue(Aallele, Aallele);
							myPopulationData.alleleC[allelePtr] = GenotypeTableUtils.getDiploidValue(Callele, Callele);
						}
					} else {
						
						//set both parents to missing if one is heterozygous
						myPopulationData.alleleA[allelePtr] = NN;
						myPopulationData.alleleC[allelePtr] = NN;
					}
				} else { //assume p0alleles.size is 2
					if (p1alleles.size() == 0) {
						myPopulationData.alleleA[allelePtr] = NN;
						myPopulationData.alleleC[allelePtr] = NN;
					} else if (p1alleles.size() == 1) {
						
						//set both parents to missing if one is heterozygous
						myPopulationData.alleleA[allelePtr] = NN;
						myPopulationData.alleleC[allelePtr] = NN;

					} else {
						myPopulationData.alleleA[allelePtr] = NN;
						myPopulationData.alleleC[allelePtr] = NN;
					}
				}
				
			}
		}
	}
	
	public ArrayList<Haplotype> mergeMajorHaplotypes(HaplotypeClusterer clusterMaker, int minClusterSize) {
		double maxMaf = 0.05;
		int maxMinorCount = 1;
		int maxDistance = 4;
		ArrayList<Haplotype> myHaplotypes = new ArrayList<Haplotype>();
		int clusterCount = 1;
		while (clusterMaker.getClusterList().get(clusterCount).getSize() >= minClusterSize) {
			Haplotype compHap = new Haplotype(clusterMaker.getClusterList().get(clusterCount).getCensoredMajorityHaplotype(maxMaf, maxMinorCount));
			for (int i = 0; i < clusterCount; i++) {
				Haplotype headHap = new Haplotype(clusterMaker.getClusterList().get(i).getCensoredMajorityHaplotype(maxMaf, maxMinorCount));
				if (headHap.distanceFrom(compHap) <= maxDistance) {
					clusterMaker.getClusterList().get(i).addAll(clusterMaker.getClusterList().get(clusterCount));
					clusterMaker.getClusterList().remove(clusterCount);
					clusterCount--;
					break;
				}
			}
			clusterCount++;
		}
		for (int i = 0; i < clusterCount; i++) myHaplotypes.add( new Haplotype(clusterMaker.getClusterList().get(i).getCensoredMajorityHaplotype(maxMaf, maxMinorCount)) );
		return myHaplotypes;
	}
	
	public void convertGenotypesToParentCalls() {
		int ntaxa = myPopulationData.original.numberOfTaxa();
		int nsites = myPopulationData.original.numberOfSites();
		byte[] parentA = myPopulationData.alleleA;
		byte[] parentC = myPopulationData.alleleC;
		GenotypeTableBuilder genoBuilder = GenotypeTableBuilder.getTaxaIncremental(myPopulationData.original.positions());
		
		for (int t = 0; t < ntaxa; t++) {
			byte[] taxongeno = myPopulationData.original.genotypeAllSites(t);
			for (int s = 0; s < nsites; s++) {
				if (taxongeno[s] != NN) {
					if (taxongeno[s] == parentA[s]) taxongeno[s] = AA;
					else if (taxongeno[s] == parentC[s]) taxongeno[s] = CC;
					else if (GenotypeTableUtils.isHeterozygous(taxongeno[s]) && parentA[s] != NN &&parentC[s] != NN) taxongeno[s] = AC;
					else taxongeno[s] = NN;
				}
			}
			genoBuilder.addTaxon(myPopulationData.original.taxa().get(t), taxongeno);
		}
		
		myPopulationData.imputed = genoBuilder.build();
	}

	public GenotypeTable preFilterSites() {
		int ntaxa = myPopulationData.original.numberOfTaxa();
//		int minCount = (int) Math.floor(ntaxa * minCoverage);
//		GenotypeTable filterGeno = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(myPopulationData.original, minMaf, 1, minCount);
		GenotypeTable filterGeno = NucleotideImputationUtils.filterSnpsByTag(myPopulationData.original, minMaf, 1 - minCoverage, 1);
		
		//select sites with proportion of hets < mean proportion het + 3 * sd(proportion het)
		int nsites = filterGeno.numberOfSites();
		double[] pHet = new double[nsites];

		for (int s = 0; s < nsites; s++) {
			int notmiss = filterGeno.totalNonMissingForSite(s);
			int hetcount = filterGeno.heterozygousCount(s);
			pHet[s] = ((double) hetcount) / notmiss;
		}

		double meanPhet = StatUtils.mean(pHet);
		double sdPhet = Math.sqrt(StatUtils.variance(pHet));
		double maxPhet = meanPhet + maxHetDeviation * sdPhet;
		
		boolean[] selected = new boolean[nsites];
		int numberFalse = 0;
		for (int s = 0; s < nsites; s++) {
			int notmiss = filterGeno.totalNonMissingForSite(s);
			int hetcount = filterGeno.heterozygousCount(s);
			double fractionHet = ((double) hetcount) / notmiss;
			if (fractionHet <= maxPhet) selected[s] = true;
			else {
				selected[s] = false;
				numberFalse++;
			}
		}
		
		//select bi-allelic sites
		int numberNotBiallelic = 0;
		for (int s = 0; s < nsites; s++) {
			int nalleles = filterGeno.alleles(s).length;
			if (nalleles != 2) {
				selected[s] = false;
				numberNotBiallelic++;
			}
		}

		
		//select sites based on minR2
		if (minR2 > 0) {
			int ldwindow = 50;
			LinkageDisequilibrium myLD = new LinkageDisequilibrium(filterGeno, ldwindow, LinkageDisequilibrium.testDesign.SlidingWindow, -1, null, false, 0, null, LinkageDisequilibrium.HetTreatment.Homozygous);
			myLD.run();
			for (int s = 0; s < nsites; s++) {
				if (selected[s]) {
					int start = Math.max(0, s - window);
					int end = Math.min(nsites, s + window);
					double sum = 0;
					double count = 0;
					for (int i = start; i <= end; i++) {
						double r2 = myLD.getRSqr(s, i);
						if (!Double.isNaN(r2)) {
							sum += myLD.getRSqr(s, i);
							count++;
						}
					}
					double r2 = sum / count;
					if (r2 < minR2) selected[s] = false;
					
				}
			}
		}
		
		int[] selectedSites = new int[nsites];
		int nSelectedSites = 0;
		for (int s = 0; s < nsites; s++) {
			if (selected[s]) selectedSites[nSelectedSites++] = s;
		}
		selectedSites = Arrays.copyOf(selectedSites, nSelectedSites);
		System.out.printf("%d sites in filtered genotype set.\n", nSelectedSites);
		return FilterGenotypeTable.getInstance(filterGeno, selectedSites);
	}
	
	//static helper functions -----------------------------------
	public static ArrayList<ArrayList<Haplotype>> getParentHaplotypesV1(ArrayList<ArrayList<Haplotype>> previousParents, ArrayList<Haplotype> candidates, int overlap, boolean forward) {
		ArrayList<ArrayList<Haplotype>> parentHaplotypes = new ArrayList<ArrayList<Haplotype>>(2);
		parentHaplotypes.add(new ArrayList<Haplotype>());
		parentHaplotypes.add(new ArrayList<Haplotype>());
		
		for (Haplotype hap : candidates) {
			boolean match = false;
			for (Haplotype parentHap : previousParents.get(0)) {
				if (doesOverlapMatch(parentHap, hap, overlap, forward)) {
					match = true;
					parentHaplotypes.get(0).add(hap);
					break;
				}
			}
			if (!match) {
				for (Haplotype parentHap : previousParents.get(1)) {
					if (doesOverlapMatch(parentHap, hap, overlap, forward)) {
						match = true;
						parentHaplotypes.get(1).add(hap);
						break;
					}
				}
			}
			if (!match) {
				StringBuilder msgBuilder = new StringBuilder("No parent matching haplotype ");
				msgBuilder.append(hap).append(" for parents:\n");
				for (Haplotype h:previousParents.get(0)) msgBuilder.append(h).append("\n");
				for (Haplotype h:previousParents.get(1)) msgBuilder.append(h).append("\n");
				System.out.println(msgBuilder.toString());
//				throw new RuntimeException(msgBuilder.toString());
				
				//find the best match
				int haplen = hap.seqlen;
				Haplotype hapstart = new Haplotype(Arrays.copyOf(hap.seq, overlap));
				int par0dist = 1000000;
				for (Haplotype parentHap : previousParents.get(0)) {
					Haplotype parentEnd = new Haplotype(Arrays.copyOfRange(parentHap.seq, haplen - overlap, haplen));
					int distFromParent = hapstart.distanceFrom(parentEnd);
					par0dist = Math.min(par0dist, distFromParent);
				}
				int par1dist = 1000000;
				for (Haplotype parentHap : previousParents.get(1)) {
					Haplotype parentEnd = new Haplotype(Arrays.copyOfRange(parentHap.seq, haplen - overlap, haplen));
					int distFromParent = hapstart.distanceFrom(parentEnd);
					par1dist = Math.min(par1dist, distFromParent);
				}
				
				if (par0dist < par1dist) parentHaplotypes.get(0).add(hap);
				else if (par0dist > par1dist) parentHaplotypes.get(1).add(hap);
				else {
					System.out.println("Haplotype not added because equi-distant from parents");
				}
			}
		}

		return parentHaplotypes;
	}
	
	public static ArrayList<ArrayList<Haplotype>> getParentHaplotypes(ArrayList<ArrayList<Haplotype>> previousParents, ArrayList<Haplotype> candidates, int overlap, boolean forward) {
		ArrayList<ArrayList<Haplotype>> parentHaplotypes = new ArrayList<ArrayList<Haplotype>>(2);
		parentHaplotypes.add(new ArrayList<Haplotype>());
		parentHaplotypes.add(new ArrayList<Haplotype>());
		int[] nPrev = new int[]{previousParents.get(0).size(), previousParents.get(1).size()};
		int maxNumber = Math.max(nPrev[0], nPrev[1]);
		
		for (Haplotype hap : candidates) {
			boolean match = false;
			for (int i = 0; i < maxNumber; i++) {
				if (i < nPrev[0] && doesOverlapMatch(previousParents.get(0).get(i), hap, overlap, forward)) {
					match = true;
					parentHaplotypes.get(0).add(hap);
					break;
				}
				if (i < nPrev[1] && doesOverlapMatch(previousParents.get(1).get(i), hap, overlap, forward)) {
					match = true;
					parentHaplotypes.get(1).add(hap);
					break;
				}
			}
			if (!match) {
				StringBuilder msgBuilder = new StringBuilder("No parent matching haplotype ");
				msgBuilder.append(hap).append(" for parents:\n");
				for (Haplotype h:previousParents.get(0)) msgBuilder.append(h).append("\n");
				for (Haplotype h:previousParents.get(1)) msgBuilder.append(h).append("\n");
				System.out.println(msgBuilder.toString());
//				throw new RuntimeException(msgBuilder.toString());
				
				//find the best match
				int haplen = hap.seqlen;
				Haplotype hapstart = new Haplotype(Arrays.copyOf(hap.seq, overlap));
				int par0dist = 1000000;
				for (Haplotype parentHap : previousParents.get(0)) {
					Haplotype parentEnd = new Haplotype(Arrays.copyOfRange(parentHap.seq, haplen - overlap, haplen));
					int distFromParent = hapstart.distanceFrom(parentEnd);
					par0dist = Math.min(par0dist, distFromParent);
				}
				int par1dist = 1000000;
				for (Haplotype parentHap : previousParents.get(1)) {
					Haplotype parentEnd = new Haplotype(Arrays.copyOfRange(parentHap.seq, haplen - overlap, haplen));
					int distFromParent = hapstart.distanceFrom(parentEnd);
					par1dist = Math.min(par1dist, distFromParent);
				}
				
				if (par0dist < par1dist) parentHaplotypes.get(0).add(hap);
				else if (par0dist > par1dist) parentHaplotypes.get(1).add(hap);
				else {
					System.out.println("Haplotype not added because equi-distant from parents");
				}
			}
		}

		return parentHaplotypes;
	}
	
	public static boolean doesOverlapMatch(Haplotype h0, Haplotype h1, int overlap, boolean forward) {
		String prev;
		if (forward) prev = h0.toString().substring(h0.seqlen - overlap, h0.seqlen);
		else prev = h0.toString().substring(0, overlap);
		
		String next;
		if (!forward) next = h1.toString().substring(h1.seqlen - overlap, h1.seqlen);
		else next = h1.toString().substring(0, overlap);
		
		boolean equals = true;
		int mismatches = 0;
		for (int i = 0; i < overlap; i++) {
			char prevChar = prev.charAt(i);
			if (prevChar != 'A' && prevChar != 'C' && prevChar != 'G' && prevChar != 'T' && prevChar != '-' && prevChar != '+') prevChar = 'N';
			char nextChar = next.charAt(i);
			if (nextChar != 'A' && nextChar != 'C' && nextChar != 'G' && nextChar != 'T' && nextChar != '-' && nextChar != '+') nextChar = 'N';
			
			if (prev.charAt(i) != next.charAt(i) && prevChar != 'N' && nextChar != 'N') {
				equals = false;
				mismatches++;
			}
		}
		return (mismatches < 2);
	}

}
