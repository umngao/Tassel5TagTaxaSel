package net.maizegenetics.analysis.imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.BinomialDistribution;

import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.Taxon;


public class RephaseParents {
	private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
	private static final byte N = GenotypeTable.UNKNOWN_ALLELE;
	private static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	private static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	private static final byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("GG");
	private static final byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("TT");
	private static final byte missingState = (byte) 4;
	
	GenotypeTable origGeno;
	Map<String, byte[]> progenyStates;
	List<String[]> plotList;
//	Map<String, String[]> plotMap;
	Map<String, byte[][]> rephasedParents;
	Map<String, byte[][]> startingParents = null;
	Map<String, double[][]> parentHaplotypeProbabilities = null; //probability that the haplotype equals the major allele
	
	//parameters
	int minFamilySize = 10;
	
	public RephaseParents() {
		
	}
	
	/**
	 * @param originalGenotypes		the GenotypeTable containing the original data being imputed
	 * @param phasedProgeny			the previously imputed progeny states (0 - 4, 4 = missing)
	 * @param plotList				the parentage information
	 * @param parentHapmap			the parent haplotypes assigned in the previous round
	 */
	public RephaseParents(GenotypeTable originalGenotypes, Map<String, byte[]> phasedProgeny, List<String[]> plotList, Map<String, byte[][]> parentHapmap) {
		origGeno = originalGenotypes;
		progenyStates = phasedProgeny;
		this.plotList = plotList;
		startingParents = parentHapmap;
		
	}
	
	Map<String, byte[][]> rephaseUsingSelfProgeny() {
		rephasedParents = new HashMap<>();
//		Map<String, List<String[]>> parentToPlotMap = plotMap.values().stream()
//				.filter(p -> p[3].equals("self"))
//				.collect(Collectors.groupingBy(p -> p[1]));
		Map<String, List<String[]>> parentToPlotMap = plotList.stream()
				.filter(p -> p[3].equals("self"))
				.collect(Collectors.groupingBy(p -> p[1]));
		
		List<String> parentList = new ArrayList<>(parentToPlotMap.keySet());
		
		//for each parent with 10 or more progeny
		int nsites = origGeno.numberOfSites();
		for (String parent : parentList) {
			List<String[]> plots = parentToPlotMap.get(parent);
			if (plots.size() >= minFamilySize) {
				//for each site
				//h0h0 progeny have h0 allele (only/mostly)
				//h1h1 progeny have h1 allele
				byte[][] newhaps = new byte[2][nsites];
				Arrays.fill(newhaps[0], N);
				rephasedParents.put(parent, newhaps);
				for (int s = 0; s < nsites; s++) {
					//for selfs, 0 and 3 are homozygotes, 1 and 2 are hets 
					//0 = h0h0, 3 = h1h1
					byte major = origGeno.majorAllele(s);
					byte minor = origGeno.minorAllele(s);
					int[] h0Counts = new int[2];
					int[] h1Counts = new int[2];
					for (String[] plot : plots) {
						//get allele counts (major, minor) for each state
						
						String progenyName = plot[0];
						int progenyIndex = origGeno.taxa().indexOf(progenyName);
						byte state = progenyStates.get(progenyName)[s];
						if (state == 0) {
							byte[] alleles = origGeno.genotypeArray(progenyIndex, s);
							if (alleles[0] == major) h0Counts[0]++;
							else if (alleles[0] == minor) h0Counts[1]++;
							if (alleles[1] == major) h0Counts[0]++;
							else if (alleles[1] == minor) h0Counts[1]++;
						}
						else if (state == 3) {
							byte[] alleles = origGeno.genotypeArray(progenyIndex, s);
							if (alleles[0] == major) h1Counts[0]++;
							else if (alleles[0] == minor) h1Counts[1]++;
							if (alleles[1] == major) h1Counts[0]++;
							else if (alleles[1] == minor) h1Counts[1]++;
						}
					}
					//if a state has only 1 allele assign that to the haplotype
					//if two alleles, if the larger count > 5 * smaller count, assign larger count to haplotype
					//this criteria might be too liberal
					int multiplier = 5;
					if (h0Counts[0] > multiplier * h0Counts[1]) newhaps[0][s] = major;
					else if (h0Counts[1] > multiplier * h0Counts[0]) newhaps[0][s] = minor;
					if (h1Counts[0] > multiplier * h1Counts[1]) newhaps[1][s] = major;
					else if (h1Counts[1] > multiplier * h1Counts[0]) newhaps[1][s] = minor;
				}
			}
		}
		
		return rephasedParents;
	}
	
	Map<String, byte[][]> rephaseUsingCrossProgenyV1() {
		int minDepth = 8;
		rephasedParents = new HashMap<>();
		
		Set<String> parentSet = plotList.stream().map(p -> p[1]).collect(Collectors.toCollection(TreeSet::new));
		parentSet.addAll(plotList.stream().map(p -> p[2]).collect(Collectors.toCollection(TreeSet::new)));
		
		Map<String, List<String[]>> parentPlotMap = new HashMap<>();
		for (String[] plot :plotList) {
			if (plot[3].equals("outcross")) {
				List<String[]> plotList = parentPlotMap.get(plot[1]);
				if (plotList == null) {
					plotList = new ArrayList<>();
					parentPlotMap.put(plot[1], plotList);
				}
				plotList.add(plot);
				if (!plot[2].equals(plot[1])) {
					plotList = parentPlotMap.get(plot[2]);
					if (plotList == null) {
						plotList = new ArrayList<>();
						parentPlotMap.put(plot[2], plotList);
					}
					plotList.add(plot);
				}
			}
		}
		
		int nsites = origGeno.numberOfSites();
		for (String parent : parentSet) {
			byte[][] newhaps = new byte[2][nsites];
			Arrays.fill(newhaps[0], N);
			Arrays.fill(newhaps[1], N);
			rephasedParents.put(parent, newhaps);
			List<byte[][]> haplotypeList = new ArrayList<>();
			List<String[]> plotList = parentPlotMap.get(parent);
			int parentIndex = origGeno.taxa().indexOf(parent);
			for (String[] plot : plotList) {
				String otherParent;
				if (plot[1].equals(parent)) otherParent = plot[2];
				else otherParent = plot[1];
				
				int progenyIndex = origGeno.taxa().indexOf(plot[0]);
				int otherIndex = origGeno.taxa().indexOf(otherParent);
				byte[] parentHap = new byte[nsites];
				byte[] myStates = progenyStates.get(plot[0]);
				for (int s = 0;s < nsites; s++) {
					byte parentGeno = origGeno.genotype(parentIndex, s);
					if (GenotypeTableUtils.isHomozygous(parentGeno) && origGeno.depth().depth(parentIndex, s) >= minDepth) {
						parentHap[s] = GenotypeTableUtils.getDiploidValues(parentGeno)[0];
					} else {
						byte progenyGeno = origGeno.genotype(progenyIndex, s);
						if (GenotypeTableUtils.isHomozygous(progenyGeno) && origGeno.depth().depth(progenyIndex, s) >= minDepth) {
							parentHap[s] = GenotypeTableUtils.getDiploidValues(progenyGeno)[0];
						} else {
							byte otherGeno = origGeno.genotype(otherIndex, s);  //maybe use other parent haplotype here instead of genotype
							if (GenotypeTableUtils.isHomozygous(otherGeno) && origGeno.depth().depth(otherIndex, s) >= minDepth && progenyGeno != NN) {
								byte[] progenyAlleles = GenotypeTableUtils.getDiploidValues(progenyGeno);
								byte[] otherAlleles = GenotypeTableUtils.getDiploidValues(otherGeno);
								if (progenyAlleles[0] != otherAlleles[0]) parentHap[s] = progenyAlleles[0];
								else if (progenyAlleles[1] != otherAlleles[0]) parentHap[s] = progenyAlleles[1];
							} else {
								//use other parent haplotype here 
								if (myStates[s] < 4 && progenyGeno != NN) {
									byte[][] otherHaps = startingParents.get(otherParent);
									//TODO finish this or remove it
								} else parentHap[s] = N;
							}
						}

					} 
				}
				haplotypeList.add(new byte[][]{parentHap, myStates});
			}
			
			//merge the haplotypes if all the same or all but one is the same(?)
			//have to use phase info here to decide which haplotypes to merge at each site to get the two parent haplotypes
			//haplotype list contains a series of inferred parent haplotypes - which haplotype? depends on site
			int nhaps = haplotypeList.size();
			int mult = 5;
			for (int s = 0;s < nsites; s++) {
				int[][] alleleCounts = new int[2][2]; //first dimension is haplotype (0 or 1), second dimension is allele (major or minor)
				byte major = origGeno.majorAllele(s);
				byte minor = origGeno.minorAllele(s);
				for (byte[][] hapinfo : haplotypeList) {
					byte state = hapinfo[1][s];
					byte hapval = hapinfo[0][s];
					if (state < 4 && hapval != N) {
						if (state < 2) {
							if (hapval == major) alleleCounts[0][0]++;
							else if (hapval == minor) alleleCounts[0][1]++;
						} else { 
							if (hapval == major) alleleCounts[1][0]++;
							else if (hapval == minor) alleleCounts[1][1]++;
						}
					}
				}
				
				//if each haplotype is mostly major or minor assign it that value
				for (int i = 0; i < 2; i++) {
					if (alleleCounts[i][0] > mult * alleleCounts[i][1]) newhaps[i][s] = major;
					else if (alleleCounts[i][1] > mult * alleleCounts[i][0]) newhaps[i][s] = minor;
				}
				
			}
			
		}
		
		return rephasedParents;
	}
	
	Map<String, byte[][]> rephaseUsingCrossProgeny() {
		int minDepth = 8;
		rephasedParents = new HashMap<>();
		Set<String> parentSet = plotList.stream().map(p -> p[1]).collect(Collectors.toCollection(TreeSet::new));
		parentSet.addAll(plotList.stream().map(p -> p[2]).collect(Collectors.toCollection(TreeSet::new)));
		
		Map<String, List<String[]>> parentPlotMap = new HashMap<>();
		for (String[] plot : plotList) {
			if (plot[3].equals("outcross")) {
				List<String[]> plotList = parentPlotMap.get(plot[1]);
				if (plotList == null) {
					plotList = new ArrayList<>();
					parentPlotMap.put(plot[1], plotList);
				}
				plotList.add(plot);
				if (!plot[2].equals(plot[1])) {
					plotList = parentPlotMap.get(plot[2]);
					if (plotList == null) {
						plotList = new ArrayList<>();
						parentPlotMap.put(plot[2], plotList);
					}
					plotList.add(plot);
				}
			}
		}
		
		int nsites = origGeno.numberOfSites();
		for (String parent : parentSet) {
			byte[][] newhaps = new byte[2][nsites];
			Arrays.fill(newhaps[0], N);
			Arrays.fill(newhaps[1], N);
			rephasedParents.put(parent, newhaps);
			List<byte[][]> haplotypeList = new ArrayList<>();
			List<String[]> plotList = parentPlotMap.get(parent);
			int parentIndex = origGeno.taxa().indexOf(parent);
			for (String[] plot : plotList) {
				String otherParent;
				boolean isFirstParent;
				if (plot[1].equals(parent)) {
					otherParent = plot[2];
					isFirstParent = true;
				}
				else {
					otherParent = plot[1];
					isFirstParent = false;
				}
				
				int progenyIndex = origGeno.taxa().indexOf(plot[0]);
				int otherIndex = origGeno.taxa().indexOf(otherParent);
				byte[][] parentHap = new byte[2][nsites];
				for (int i = 0; i < 2; i++) Arrays.fill(parentHap[i], N);
				byte[] myStates = progenyStates.get(plot[0]);
				
				for (int s = 0;s < nsites; s++) {
					byte siteState = myStates[s];
					if (siteState == 4) continue; //unknown state, do not evaluate this site
					
					//at this site, which haplotype did this parent pass on to this progeny
					int progenyHaplotypeIndex, nonProgenyHaplotypeIndex;
					if (isFirstParent) {
						if (myStates[s] == 0 || myStates[s] == 1) {
							progenyHaplotypeIndex = 0;
							nonProgenyHaplotypeIndex = 1;
						} else {
							progenyHaplotypeIndex = 1;
							nonProgenyHaplotypeIndex = 0;
						}
					} else {
						if (myStates[s] == 0 || myStates[s] == 2) {
							progenyHaplotypeIndex = 0;
							nonProgenyHaplotypeIndex = 1;
						} else {
							progenyHaplotypeIndex = 1;
							nonProgenyHaplotypeIndex = 0;
						}
					}
					
					//test homozygosity by combination of isHomozygous and minimum depth
					byte parentGeno = origGeno.genotype(parentIndex, s);
					if (GenotypeTableUtils.isHomozygous(parentGeno) && origGeno.depth().depth(parentIndex, s) >= minDepth) {
						//if the parent is homozygous then it has the same allele in both haplotypes
						parentHap[0][s] = parentHap[1][s] = GenotypeTableUtils.getDiploidValues(parentGeno)[0];
					} else {
						byte progenyGeno = origGeno.genotype(progenyIndex, s);
						if (GenotypeTableUtils.isHomozygous(progenyGeno) && origGeno.depth().depth(progenyIndex, s) >= minDepth) {
							//if progeny is homozygous then the parent haplotype that it carries has that allele
							parentHap[progenyHaplotypeIndex][s] = GenotypeTableUtils.getDiploidValues(progenyGeno)[0];
							if (parentGeno != NN ) {
								//if the parent has an allele other than the progeny allele, than its non-progeny haplotype is the non-progeny allele
								byte[] progenyAlleles = GenotypeTableUtils.getDiploidValues(progenyGeno);
								byte[] parentAlleles = GenotypeTableUtils.getDiploidValues(parentGeno);
								if (progenyAlleles[0] != parentAlleles[0]) parentHap[nonProgenyHaplotypeIndex][s] = progenyAlleles[0];
								else if (progenyAlleles[1] != parentAlleles[0]) parentHap[nonProgenyHaplotypeIndex][s] = progenyAlleles[1];
							}
						} else {
							byte otherGeno = origGeno.genotype(otherIndex, s);
							//if other Parent is homozygous and progeny carries the alternate allele then the parent must have contributed that
							if (GenotypeTableUtils.isHomozygous(otherGeno) && origGeno.depth().depth(otherIndex, s) >= minDepth && progenyGeno != NN) {
								byte[] progenyAlleles = GenotypeTableUtils.getDiploidValues(progenyGeno);
								byte[] otherAlleles = GenotypeTableUtils.getDiploidValues(otherGeno);
								if (progenyAlleles[0] != otherAlleles[0] || progenyAlleles[1] != otherAlleles[0]) {
									if (progenyAlleles[0] != otherAlleles[0]) parentHap[progenyHaplotypeIndex][s] = progenyAlleles[0];
									else parentHap[progenyHaplotypeIndex][s] = progenyAlleles[1];
									//in addition, if the parent carries the otherParent allele then the other haplotype must have that allele
									if (parentGeno == otherGeno || GenotypeTableUtils.isHeterozygous(parentGeno)) parentHap[nonProgenyHaplotypeIndex][s] = otherAlleles[0];
								}
							} else {
								//if other parent is not known to be homozygous, use the other parent haplotype if known
								if (myStates[s] < 4 && progenyGeno != NN) {
									byte[][] otherHaps = startingParents.get(otherParent);
									byte otherParentHaplotype;
									if (progenyHaplotypeIndex == 0) {
										if (myStates[s] == 0 || myStates[s] == 2) otherParentHaplotype = otherHaps[0][s];
										else otherParentHaplotype = otherHaps[1][s];
									} else {
										if (myStates[s] == 0 || myStates[s] == 1) otherParentHaplotype = otherHaps[0][s];
										else otherParentHaplotype = otherHaps[1][s];
									}
									
									if (otherParentHaplotype != N) {
										//repeat logic for other parent homozygous
										byte[] progenyAlleles = GenotypeTableUtils.getDiploidValues(progenyGeno);
										if (progenyAlleles[0] != otherParentHaplotype || progenyAlleles[1] != otherParentHaplotype) {
											if (progenyAlleles[0] != otherParentHaplotype) parentHap[progenyHaplotypeIndex][s] = progenyAlleles[0];
											else parentHap[progenyHaplotypeIndex][s] = progenyAlleles[1];
											byte[] parentAlleles = GenotypeTableUtils.getDiploidValues(parentGeno);
											if (parentAlleles[0] != otherParentHaplotype || parentAlleles[1] != otherParentHaplotype) parentHap[nonProgenyHaplotypeIndex][s] = otherParentHaplotype;
										}
									}
								}
								
								
							}
						}

					} 
				}
				haplotypeList.add(parentHap);
			}
			
			//merge the haplotypes if all the same or all but one is the same(?)
			//have to use phase info here to decide which haplotypes to merge at each site to get the two parent haplotypes
			//haplotype list contains a series of inferred parent haplotypes - which haplotype? depends on site
			int nhaps = haplotypeList.size();
			int mult = 5;
			for (int s = 0;s < nsites; s++) {
				int[][] alleleCounts = new int[2][2]; //first dimension is haplotype (0 or 1), second dimension is allele (major or minor)
				byte major = origGeno.majorAllele(s);
				byte minor = origGeno.minorAllele(s);
				for (byte[][] hapinfo : haplotypeList) {
					byte state = hapinfo[1][s];
					byte hapval = hapinfo[0][s];
					if (state < 4 && hapval != N) {
						if (state < 2) {
							if (hapval == major) alleleCounts[0][0]++;
							else if (hapval == minor) alleleCounts[0][1]++;
						} else { 
							if (hapval == major) alleleCounts[1][0]++;
							else if (hapval == minor) alleleCounts[1][1]++;
						}
					}
				}
				
				//if each haplotype is mostly major or minor assign it that value
				for (int i = 0; i < 2; i++) {
					if (alleleCounts[i][0] > mult * alleleCounts[i][1]) newhaps[i][s] = major;
					else if (alleleCounts[i][1] > mult * alleleCounts[i][0]) newhaps[i][s] = minor;
				}
				
			}
			
		}
		
		return rephasedParents;
	}
	
	Map<String, double[][]> rephaseUsingAlleleDepth() {
		return rephaseUsingAlleleDepth(null);
	}
	
	Map<String, double[][]> rephaseUsingAlleleDepth(String saveFilename) {
		Map<String, double[][]> nextHaplotypeProbs = new HashMap<>();

		//can phase each parent using all progeny with a phased other parent
		//for each parent create a list of all plots for that parent for which the other parent is phased
		Map<String, List<String[]>> parentPlotMap = new HashMap<>();
		for (String[] plot : plotList) {

			if (startingParents.get(plot[2]) != null) {
				List<String[]> parentPlotList = parentPlotMap.get(plot[1]);
				if (parentPlotList == null) {
					parentPlotList = new ArrayList<>();
					parentPlotMap.put(plot[1], parentPlotList);
				}
				parentPlotList.add(plot);
			}
			
			if (!plot[2].equals(plot[1]) && startingParents.get(plot[1]) != null) {
				List<String[]> parentPlotList = parentPlotMap.get(plot[2]);
				if (parentPlotList == null) {
					parentPlotList = new ArrayList<>();
					parentPlotMap.put(plot[2], parentPlotList);
				}
				parentPlotList.add(plot);
			}

		}
		
		//for each parent haplotype, at each site, calculate the probability that the haplotype is the major allele
		for (String parent : parentPlotMap.keySet()) {
			if (startingParents.get(parent) != null) {
				double[][] probs = rephasePreviouslyPhased(parent, parentPlotMap.get(parent));
				nextHaplotypeProbs.put(parent, probs);
			}
		}
		
		parentHaplotypeProbabilities = nextHaplotypeProbs;
		
		//save the rephased haplotypes to a file
//		String savename = "/Users/pbradbury/Documents/projects/teosinte/haplotypes/C2_rephased_parents.bin";
		if (saveFilename != null && saveFilename.length() > 1) {
			try {
				FileOutputStream fos = new FileOutputStream(new File(saveFilename));
				ObjectOutputStream oos = new ObjectOutputStream(fos);
				oos.writeObject(parentHaplotypeProbabilities);
				oos.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return parentHaplotypeProbabilities;
	}
	
	double[][] rephasePreviouslyPhasedV2(String parent, List<String[]> plotList) {
		
		double err = 0.01;
		int nsites = origGeno.numberOfSites();
		double[][] haplotypeProbability = new double[2][nsites];
		for (int i = 0; i < 2; i++) Arrays.fill(haplotypeProbability[i], Double.NaN);
		
		for (int s = 0; s < nsites; s++) {
			byte major = origGeno.majorAllele(s);
			byte minor = origGeno.minorAllele(s);
			if (minor == N) {
				for (int i = 0; i < 2; i++) haplotypeProbability[i][s] = 1;
			} else {
				double majorFreq = origGeno.majorAlleleFrequency(s);
				double minorFreq = origGeno.minorAlleleFrequency(s);

				//sum allele depths across all plots and by progeny state
				int[] totalAlleleDepths = new int[6];
				int[][] stateAlleleDepths = new int[4][6];
				int[][] parentStateAlleleDepths = new int[2][6];
				int[] switchState = new int[]{0,2,1,3};
				for (String[] plot : plotList) {
					int taxonNdx = origGeno.taxa().indexOf(plot[0]);
					int[] tempDepths = origGeno.depthForAlleles(taxonNdx, s);
					int state = progenyStates.get(plot[0])[s];
					if (state < 4) {
						int parentState;
						int myState;
						if (parent.equals(plot[1])) {
							if (state == 0 || state == 1) parentState = 0;
							else parentState = 1;
							myState = state;
						} else {
							if (state == 0 || state == 2) parentState = 0;
							else parentState = 1;
							myState = switchState[state];
						}
						for (int i = 0; i < 6; i++) {
							totalAlleleDepths[i] += tempDepths[i];
							stateAlleleDepths[myState][i] += tempDepths[i];
							parentStateAlleleDepths[parentState][i] += tempDepths[i];
						}

					}
				}
				
				if (totalAlleleDepths[major] < 10) {
					//skip, too little data
				} else {
					//P(hap=major|obs) = P(obs|hap=major)P(hap=major)/P(obs)
					//P(hap=minor|obs) = P(obs|hap=minor)P(hap=minor)/P(obs)
					//P(hap=major|obs) + P(hap=minor|obs) = 1
					//P(hap=major|obs) = P(obs|hap=major)P(hap=major) / [P(obs|hap=major)P(hap=major) + P(obs|hap=minor)P(hap=minor)/P(obs)]
					
					//P(obs|hap=major) = binom(trials = total depth, successes = major depth, (1 + maf)/2)
					//P(obs|hap=minor) = binom(total depth, minor depth, (1 + minorFreq)/2)
					
					for (int chr = 0; chr < 2; chr++) {
						int majorDepth = parentStateAlleleDepths[chr][major];
						int minorDepth = parentStateAlleleDepths[chr][minor];
						int totalDepth = majorDepth + minorDepth;
						double pmajor = new BinomialDistribution(totalDepth, (1 + majorFreq)/2).probability(majorDepth);
						pmajor *= majorFreq;
						double pminor = new BinomialDistribution(totalDepth, (1 + minorFreq)/2).probability(majorDepth);
						pminor *= minorFreq;
						haplotypeProbability[chr][s] = pmajor / (pmajor + pminor);
					}
					
					
				}
				
			}
		}
		return haplotypeProbability;
	}
	
	double[][] rephasePreviouslyPhased(String parent, List<String[]> plotList) {
		//calculate P(hap=major) for each haplotype
		//P(hap=a | obs) = P(hap=a | obs,other=a)*P(other=a) + P(hap=a | obs,other=b)*P(other=b)
		//	~ P(obs | hap=a, other=a)*P(other=a) + P(obs | hap=a, other=b)*P(other=b)
		
		//for each site
		//P(hap=a)=P(hap=b)=0.5
		//for P(other=a), P(other=b) = allele freq of a and b in total population
		//for each parent haplotype
		//group progeny into those with other parent carrying the major allele (a) vs. the minor allele (b)
		//add up the allele depths for each group. These are the obs.
		//calculate P(obs|...) from the binomial distribution
		//calculate P(hap=major) and store the result
		
		double err = 0.01;
		int nsites = origGeno.numberOfSites();
		double[][] haplotypeProbability = new double[2][nsites];
		for (int i = 0; i < 2; i++) Arrays.fill(haplotypeProbability[i], Double.NaN);
		
		for (int s = 0; s < nsites; s++) {
			byte major = origGeno.majorAllele(s);
			byte minor = origGeno.minorAllele(s);
			if (minor == N) {
				for (int i = 0; i < 2; i++) haplotypeProbability[i][s] = 1;
			} else {
				double majorFreq = origGeno.majorAlleleFrequency(s);
				double minorFreq = origGeno.minorAlleleFrequency(s);

				//sum allele depths across all plots and by progeny state
				int[] totalAlleleDepths = new int[6];
				int[][] stateAlleleDepths = new int[4][6];
				int[][] parentStateAlleleDepths = new int[2][6];
				int[] switchState = new int[]{0,2,1,3};
				for (String[] plot : plotList) {
					int taxonNdx = origGeno.taxa().indexOf(plot[0]);
					int[] tempDepths = origGeno.depthForAlleles(taxonNdx, s);
					byte geno = origGeno.genotype(taxonNdx, s);
					int state = progenyStates.get(plot[0])[s];
					if (state < 4) {
						int parentState;
						int myState;
						if (parent.equals(plot[1])) {
							if (state == 0 || state == 1) parentState = 0;
							else parentState = 1;
							myState = state;
						} else {
							if (state == 0 || state == 2) parentState = 0;
							else parentState = 1;
							myState = switchState[state];
						}
						for (int i = 0; i < 6; i++) {
							totalAlleleDepths[i] += tempDepths[i];
							stateAlleleDepths[myState][i] += tempDepths[i];
							parentStateAlleleDepths[parentState][i] += tempDepths[i];
						}

					}
				}
				
				if (parentStateAlleleDepths[0][major] < 5 || parentStateAlleleDepths[0][major] < 5) {
					//skip, too little data
				} else {
					//for each parent haplotype
					//P(hap=a | obs) = P(obs | hap=a, other=a)*P(hap=a)*P(other=a) + P(obs | hap=a, other=b)*P(hap=a)*P(other=b)
					//then divide by sum of P(hap=a|obs) + P(hap=b|obs)
					//use origGeno allele freq for p(other=a) and p(hap=a)
					//consider using something other than allele frequency for this
					//perhaps using previously computed p's
					
					//parent haplotype 0, pMajor. Only consider the observations for h00 genotypes. h01 does not need to be considered
					//P(obs | h00=major) = [P(obs | h00=major, h10=major, h11=major) * P(h10=major) *P(h11=major)
					//                   + P(obs | h00=major, h10=major, h11=minor) * P(h10=major) * P(h11=minor)
					//                   + P(obs | h00=major, h10=minor, h11=major) * P(h10=minor) * P(h11=major)
					//                   + P(obs | h00=major, h10=minor, h11=minor) * P(h10=minor) * P(h11=minor))] * P(h00=major)

					//for h00 = major
					double ph10major = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], err).probability(stateAlleleDepths[0][minor]);  //h10 = major
					double ph10minor = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], 0.5).probability(stateAlleleDepths[0][minor]);  //h10 = minor
					double ph11major = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], err).probability(stateAlleleDepths[1][minor]);  //h11 = major
					double ph11minor = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], 0.5).probability(stateAlleleDepths[1][minor]);  //h11 = minor
					double ph00major = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph00major *= majorFreq;
					
					//for h00 = minor
					ph10major = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], 0.5).probability(stateAlleleDepths[0][major]);  //h10 = major
					ph10minor = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], err).probability(stateAlleleDepths[0][major]);  //h10 = minor
					ph11major = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], 0.5).probability(stateAlleleDepths[1][major]);  //h11 = major
					ph11minor = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], err).probability(stateAlleleDepths[1][major]);  //h11 = minor
					double ph00minor = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph00minor *= minorFreq;
					
					haplotypeProbability[0][s] = ph00major / (ph00major + ph00minor);
					
					//for h01 = major
					ph10major = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], err).probability(stateAlleleDepths[2][minor]);  //h10 = major
					ph10minor = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], 0.5).probability(stateAlleleDepths[2][minor]);  //h10 = minor
					ph11major = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], err).probability(stateAlleleDepths[3][minor]);  //h11 = major
					ph11minor = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], 0.5).probability(stateAlleleDepths[3][minor]);  //h11 = minor
					double ph01major = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph01major *= majorFreq;
					
					//for h01 = minor
					ph10major = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], 0.5).probability(stateAlleleDepths[2][major]);  //h10 = major
					ph10minor = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], err).probability(stateAlleleDepths[2][major]);  //h10 = minor
					ph11major = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], 0.5).probability(stateAlleleDepths[3][major]);  //h11 = major
					ph11minor = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], err).probability(stateAlleleDepths[3][major]);  //h11 = minor
					double ph01minor = ph10major * ph11major * majorFreq * majorFreq 
							+ ph10major * ph11minor * majorFreq * minorFreq 
							+ ph10minor * ph11major * minorFreq * majorFreq 
							+ ph10minor * ph11minor * minorFreq * minorFreq;
					ph01minor *= minorFreq;
					
					haplotypeProbability[1][s] = ph01major / (ph01major + ph01minor);

				}
				
			}
		}
		return haplotypeProbability;
	}
	
	double[][] rephasePreviouslyPhasedMultithreaded(String parent, List<String[]> progenyPlotList) {
		int maxThreads = 10;
		int nThreads = Runtime.getRuntime().availableProcessors() - 1;
		nThreads = Math.min(nThreads, maxThreads);
		ExecutorService exec = Executors.newFixedThreadPool(nThreads);
		List<HaplotypeProbabilityThread> myThreads = new ArrayList<>();
		BlockingQueue<Optional<SitePlotValues>> inputQueue = new LinkedBlockingQueue<>(500);
		BlockingQueue<Optional<SitePlotValues>> outputQueue = new LinkedBlockingQueue<>();
		
		//create the worker threads
		for (int t = 0; t < nThreads; t++) {
			myThreads.add(new HaplotypeProbabilityThread(parent, progenyPlotList, inputQueue, outputQueue));
		}

		//start the worker threads
		for (Thread thr : myThreads) exec.submit(thr);
		exec.shutdown();
		
		//submit sites to the worker threads
		int nsites = origGeno.numberOfSites();
		String[] taxanames = progenyPlotList.stream().map(p -> p[0]).toArray(String[]::new);
		Arrays.sort(taxanames);
		int[] taxaIndex = Arrays.stream(taxanames).mapToInt(tn -> origGeno.taxa().indexOf(tn)).toArray();
		
		for (int s = 0; s < nsites; s++) {
			SitePlotValues spv = new SitePlotValues(origGeno, progenyStates, s, taxanames, taxaIndex);
			boolean successful = false;
			try {
				successful = inputQueue.offer(Optional.of(spv), 5, TimeUnit.SECONDS);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			if (!successful) {
				System.out.printf("Failed to add site %d to input queue.\n", s);
			}
		}

		for (int i = 0; i < nThreads; i++) {
			try {
				boolean successful = inputQueue.offer(Optional.empty(), 1, TimeUnit.SECONDS);
				if (!successful) {
					System.out.println("Failed to add an end signal (empty Optional) to the input queue.");
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		//process the output
		double[][] haplotypeProbability = new double[2][nsites];
		for (int i = 0; i < 2; i++) Arrays.fill(haplotypeProbability[i], Double.NaN);

		int numberOfThreadsFinished = 0;
		while (!exec.isTerminated()) {
			try {
				Optional<SitePlotValues> opt = outputQueue.poll(1, TimeUnit.SECONDS);
				if (opt != null && opt.isPresent()) {
					SitePlotValues spv = opt.get();
					int site = spv.site;
					haplotypeProbability[0][site] = spv.prob[0];
					haplotypeProbability[1][site] = spv.prob[1];
				} 
			} catch (InterruptedException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		return haplotypeProbability;
	}
	
	double[] siteHaplotypeProbabilities(String parent, List<String[]> progenyPlotList, SitePlotValues spv) {
		double err = 0.01;
		double[] siteProbs = new double[2];
		if (spv.minor == N) {
			for (int i = 0; i < 2; i++) siteProbs[i] = 1;
		} else {
			//sum allele depths across all plots and by progeny state
			//this code assumes that A = parent1 chr0 + parent2 chr0
			//C = parent1 chr0 + parent2 chr1
			//G = parent1 chr1 + parent2 chr0
			//T = parent1 chr1 + parent2 chr1
			//if parent = parent 2 then C and G have to be switched to make all the calls consistent
			int major = spv.major;
			int minor = spv.minor;
			double majorFreq = spv.majorFreq;
			double minorFreq = spv.minorFreq;
			int nprogeny = progenyPlotList.size();
			int[] totalAlleleDepths = new int[6];
			int[][] stateAlleleDepths = new int[4][6];
			int[] switchState = new int[]{0,2,1,3}; //switches 1 and 2 (C and G)
			for (int p = 0; p < nprogeny; p++) {
				String[] plot = progenyPlotList.get(p);
				int state = spv.plotStates[p];
				int[] depths = spv.depthList.get(p);
				if (state < 4) {
					int myState;
					if (parent.equals(plot[1])) {
						myState = state;
					} else {
						myState = switchState[state];
					}
					for (int i = 0; i < 6; i++) {
						totalAlleleDepths[i] += depths[i];
						stateAlleleDepths[myState][i] += depths[i];
					}

				}
			}
			
			if (totalAlleleDepths[major] < 10) {
				//skip, too little data
			} else {
				//for each parent haplotype
				//P(hap=a | obs) = P(obs | hap=a, other=a)*P(hap=a)*P(other=a) + P(obs | hap=a, other=b)*P(hap=a)*P(other=b)
				//then divide by sum of P(hap=a|obs) + P(hap=b|obs)
				//use P(hap=a) = 0.5. Since a constant can drop from equation
				//use origGeno allele freq for p(other=a)
				
				//parent haplotype 0, pMajor. Only consider the observations for h00 genotypes. h01 does not need to be considered
				//P(obs | h00=major) = P(obs | h00=major, h10=major, h11=major) * P(h10=major) *P(h11=major)
				//                   + P(obs | h00=major, h10=major, h11=minor) * P(h10=major) * P(h11=minor)
				//                   + P(obs | h00=major, h10=minor, h11=major) * P(h10=minor) * P(h11=major)
				//                   + P(obs | h00=major, h10=minor, h11=minor) * P(h10=minor) * P(h11=minor)

				//for h00 = major
				double ph10major = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], err).probability(stateAlleleDepths[0][minor]);  //h10 = major
				double ph10minor = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], 0.5).probability(stateAlleleDepths[0][minor]);  //h10 = minor
				double ph11major = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], err).probability(stateAlleleDepths[1][minor]);  //h11 = major
				double ph11minor = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], 0.5).probability(stateAlleleDepths[1][minor]);  //h11 = minor
				double ph00major = ph10major * ph11major * majorFreq * majorFreq 
						+ ph10major * ph11minor * majorFreq * minorFreq 
						+ ph10minor * ph11major * minorFreq * majorFreq 
						+ ph10minor * ph11minor * minorFreq * minorFreq;
				ph00major *= majorFreq;
				
				//for h00 = minor
				ph10major = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], 0.5).probability(stateAlleleDepths[0][major]);  //h10 = major
				ph10minor = new BinomialDistribution(stateAlleleDepths[0][major] + stateAlleleDepths[0][minor], err).probability(stateAlleleDepths[0][major]);  //h10 = minor
				ph11major = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], 0.5).probability(stateAlleleDepths[1][major]);  //h11 = major
				ph11minor = new BinomialDistribution(stateAlleleDepths[1][major] + stateAlleleDepths[1][minor], err).probability(stateAlleleDepths[1][major]);  //h11 = minor
				double ph00minor = ph10major * ph11major * majorFreq * majorFreq 
						+ ph10major * ph11minor * majorFreq * minorFreq 
						+ ph10minor * ph11major * minorFreq * majorFreq 
						+ ph10minor * ph11minor * minorFreq * minorFreq;
				ph00minor *= minorFreq;
				
				siteProbs[0] = ph00major / (ph00major + ph00minor);
				
				//for h01 = major
				ph10major = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], err).probability(stateAlleleDepths[2][minor]);  //h10 = major
				ph10minor = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], 0.5).probability(stateAlleleDepths[2][minor]);  //h10 = minor
				ph11major = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], err).probability(stateAlleleDepths[3][minor]);  //h11 = major
				ph11minor = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], 0.5).probability(stateAlleleDepths[3][minor]);  //h11 = minor
				double ph01major = ph10major * ph11major * majorFreq * majorFreq 
						+ ph10major * ph11minor * majorFreq * minorFreq 
						+ ph10minor * ph11major * minorFreq * majorFreq 
						+ ph10minor * ph11minor * minorFreq * minorFreq;
				ph01major *= majorFreq;
				
				//for h01 = minor
				ph10major = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], 0.5).probability(stateAlleleDepths[2][major]);  //h10 = major
				ph10minor = new BinomialDistribution(stateAlleleDepths[2][major] + stateAlleleDepths[2][minor], err).probability(stateAlleleDepths[2][major]);  //h10 = minor
				ph11major = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], 0.5).probability(stateAlleleDepths[3][major]);  //h11 = major
				ph11minor = new BinomialDistribution(stateAlleleDepths[3][major] + stateAlleleDepths[3][minor], err).probability(stateAlleleDepths[3][major]);  //h11 = minor
				double ph01minor = ph10major * ph11major * majorFreq * majorFreq 
						+ ph10major * ph11minor * majorFreq * minorFreq 
						+ ph10minor * ph11major * minorFreq * majorFreq 
						+ ph10minor * ph11minor * minorFreq * minorFreq;
				ph01minor *= minorFreq;
				
				siteProbs[1] = ph01major / (ph01major + ph01minor);

			}
			
		}
		
		return siteProbs;
	}
	
	Map<String, byte[]> readProgenyStates() {
		Map<String, byte[]> states = new HashMap<>();
		FileLoadPlugin flp = new FileLoadPlugin(null, false);
		flp.setTheFileType(TasselFileType.Unknown);
		flp.setOpenFiles(new File[]{new File("")});
		GenotypeTable stateGeno = (GenotypeTable) flp.performFunction(null).getData(0).getData();

		
		return states;
	}

	private class SitePlotValues {
		int site;
		byte major, minor;
		double majorFreq, minorFreq;
		List<int[]> depthList;
		int[] plotStates;
		double[] prob;
		
		SitePlotValues (int site, byte major, byte minor, double majorFreq, double minorFreq, List<int[]> depthList, int[] plotStates) {
			this.site = site;
			this.major = major;
			this.minor = minor;
			this.majorFreq = majorFreq;
			this.minorFreq = minorFreq;
			this.depthList = depthList;
			this.plotStates = plotStates;
		}
		
		SitePlotValues(GenotypeTable origGeno, Map<String, byte[]> taxaStates, int site, String[] taxanames, int[] taxaIndex) {
			this.site = site;
			major = origGeno.majorAllele(site);
			minor = origGeno.minorAllele(site);
			majorFreq = origGeno.majorAlleleFrequency(site);
			minorFreq = origGeno.minorAlleleFrequency(site);
			
			int ntaxa = taxanames.length;
			depthList = new ArrayList<>();
			plotStates = new int[ntaxa];
			for(int t = 0; t < ntaxa; t++) {
				depthList.add(origGeno.depthForAlleles(taxaIndex[t], site));
				plotStates[t] = (int) taxaStates.get(taxanames[t])[site];
			}
		}
		
//		void setProbabilities(double[] probs) { prob = probs; }
	}

	private class HaplotypeProbabilityThread extends Thread {
		String parent;
		List<String[]> progenyPlotList;
		BlockingQueue<Optional<SitePlotValues>> siteQueue;
		BlockingQueue<Optional<SitePlotValues>> resultQueue;
		
		HaplotypeProbabilityThread(String parent, List<String[]> progenyPlotList, BlockingQueue<Optional<SitePlotValues>> siteQueue, BlockingQueue<Optional<SitePlotValues>> resultQueue) {
			
			this.parent = parent;
			this.progenyPlotList = progenyPlotList;
			this.siteQueue = siteQueue;
			this.resultQueue = resultQueue;
		}
		
		@Override
		public void run() {
			Optional<SitePlotValues> opt = Optional.empty();
			System.out.println("Starting new HaplotypeProbabilityThread.");
			try {
				opt = siteQueue.poll(1, TimeUnit.SECONDS);
			} catch (InterruptedException e) {
				e.printStackTrace();
				return;
			}
			
			boolean successful = true;
			while (opt != null && opt.isPresent() && successful) {
				SitePlotValues spv = opt.get();
				 
				//calculate the probability
				double[] prob = siteHaplotypeProbabilities(parent, progenyPlotList, spv);
				spv.prob = prob;
				try {
					Optional<SitePlotValues> result = Optional.of(spv);
					successful = resultQueue.offer(result, 5, TimeUnit.SECONDS);
					if (successful) {
						opt = siteQueue.poll(1, TimeUnit.SECONDS);
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
					return;
				}
				if (!successful) System.out.printf("failed to add spv to result queue at site = %d\n", spv.site);
			}
			if (opt == null) System.out.println("HaplotypeProbabilityThread finished because opt is null");
			if (!opt.isPresent()) System.out.println("HaplotypeProbabilityThread finished because opt is empty");
			System.out.println("Returning from HaplotypeProbabilityThread.");
		}
		
	}
}
