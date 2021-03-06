package net.maizegenetics.analysis.clustering;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multiset.Entry;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 * A HaplotypeCluster is a List of Haplotypes that are all similar. Size and score are
 * alternate methods of measuring cluster importance. For score, Haplotypes are weighted by the
 * number of clusters to which they belong since a Haplotype can belong to more than one cluster.
 * The natural order for clusters is determined by score then by size, both in descending order.
 * @author Peter Bradbury
 */
public class HaplotypeCluster implements Comparable<HaplotypeCluster> {
	private ArrayList<Haplotype> hapList;
	private double score = 0;
	private ArrayList<int[][]> alleleCounts = null;

	/**
	 * TYPE determines how the haplotype will be represented. 
	 * If majority, the most common allele at a locus within the cluster will be used.
	 * If unanimous, then the alleles used will be those for which almost all taxa within the cluster are the same. One error is allowed.
	 */
	public enum TYPE {majority, unanimous, censored}
	
	/**
	 * The TYPE used to represent all HaplotypeClusters
	 */
	public static TYPE ReturnHaplotype = TYPE.unanimous;
	
	/**
	 * The byte value for missing. It is set to NucleotideAlignmentConstants.getNucleotideDiploidByte('N').
	 */
	static byte N = NucleotideAlignmentConstants.getNucleotideDiploidByte('N'); 
	
	/**
	 * @param hapList	an ArrayList of Haplotypes
	 */
	public HaplotypeCluster(ArrayList<Haplotype> hapList) {
		this.hapList = hapList;
	}
	
	/**
	 * @param hap	a single Haplotype used to seed the cluster
	 */
	public HaplotypeCluster(Haplotype hap) {
		hapList = new ArrayList<Haplotype>();
		hapList.add(hap);
	}

	/**
	 * @param hap	a Haplotype
	 * @param initialScore	the initial cluster score, usually set to 1
	 */
	public HaplotypeCluster(Haplotype hap, double initialScore) {
		hapList = new ArrayList<Haplotype>();
		hapList.add(hap);
		score = initialScore;
	}

	/**
	 * @param hapList	an ArrayList of Haplotype objects
	 * @param initialScore	the initial score assigned to this cluster
	 */
	public HaplotypeCluster(ArrayList<Haplotype> hapList, double initialScore) {
		this.hapList = hapList;
		score = initialScore;
	}
	
	
	/**
	 * @return	the number of Haplotypes in this cluster
	 */
	public int getSize() { return hapList.size(); }
	
	/**
	 * @return	the Haplotypes in the cluster
	 */
	public ArrayList<Haplotype> getHaplotypeList() { return hapList; }
	
	/**
	 * @return	the cluster score. Typically this represents the importance of a cluster but be used
	 * for anything by the calling program.
	 * A typical use is as the sum of the individual Haplotype scores where a score equals
	 * the 1 divided by the number of clusters containing a Haplotype.
	 */
	public double getScore() { return score; }
	
	/**
	 * Sets the score of this haplotype to the new score.
	 * @param newScore	the new score
	 */
	public void setScore(double newScore) { score = newScore; }
	
	/**
	 * @param val	the value to be added to the cluster score
	 */
	public void incrementScore(double val) { score += val;}
	
	/**
	 * @param index	the index of a Haplotype in this cluster
	 * @return	the Haplotype indexed by index
	 */
	public Haplotype get(int index) { return hapList.get(index); }
	
	/**
	 * @param hap	a Haplotype
	 */
	public void add(Haplotype hap) { hapList.add(hap); }
	
	/**
	 * @param hap	a Haplotype to be removed from the cluster
	 */
	public void remove(Haplotype hap) { hapList.remove(hap); }
	
	/**
	 * All haplotypes in cluster that are also in this cluster will be removed from this cluster
	 * @param cluster	another cluster
	 */
	public void removeAll(HaplotypeCluster cluster) { hapList.removeAll(cluster.hapList); }
	
	/**
	 * This function does not prevent Haplotypes from being duplicated.
	 * To avoid duplicates use addAllUnique.
	 * @param cluster	Adds all Haplotypes in cluster to this cluster.
	 */
	public void addAll(HaplotypeCluster cluster) { hapList.addAll(cluster.hapList); }
	
	/**
	 * @return	an Iterator for the Haplotypes in this cluster
	 */
	public Iterator<Haplotype> getIterator() { return hapList.iterator(); }
	
	/**
	 * Adds Haplotypes from cluster that are not already contained in this cluster
	 * @param cluster	another cluster
	 */
	public void addAllUnique(HaplotypeCluster cluster) {
		Iterator<Haplotype> hit = cluster.getIterator();
		while (hit.hasNext()) {
			Haplotype hap = hit.next();
			if (!hapList.contains(hap)) hapList.add(hap);
		}
	}
	
	/**
	 * Counts the number of Haplotypes in cluster that are not also in this cluster
	 * @param cluster	another cluster
	 * @return	the count of Haplotypes unique to cluster
	 */
	public int getCountOfHaplotypesNotInThisCluster(HaplotypeCluster cluster) {
		int hapcount = 0;
		Iterator<Haplotype> hit = cluster.getIterator();
		while (hit.hasNext()) {
			Haplotype hap = hit.next();
			if (!hapList.contains(hap)) hapcount++;
		}
		return hapcount;
	}
	
	/**
	 * @param hap	a Haplotype
	 * @return the maximum distance of this Haplotype from each Haplotype in this cluster
	 */
	public int getMaximumDistance(Haplotype hap) {
		int maxd = 0;
		for (Haplotype myHaplotype : hapList) {
			int dist = myHaplotype.distanceFrom(hap);
			maxd = Math.max(maxd, dist);
		}
		return maxd;
	}
	
	/**
	 * 
	 * @param site
	 * @return	the alleles at this site with respective count, sorted by count, high to low
	 * The first dimension is the allele number. For the second dimension, the zero element is the byte value of the allele and 
	 * the second element is the number of individuals carrying that allele.
	 */
	public int[][] getAllelesAtSite(int site) {
		if (alleleCounts == null) countAllelesAtAllSites();
		return alleleCounts.get(site);
	}
	
	public void countAllelesAtAllSites() {
		alleleCounts = new ArrayList<int[][]>();
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		int nsites = hapList.get(0).seqlen;
		int nhaps = hapList.size();
		
		for (int s = 0; s < nsites; s++) {
			Multiset<Byte> alleleset = HashMultiset.create();
			for (Haplotype hap : hapList) {
				byte val = hap.seq[s];
				if (val != NN) alleleset.add(val);
			}
			int nalleles = alleleset.elementSet().size();
			int[][] alleleStats = new int[nalleles][2];
			int alleleCount = 0;
			for (Byte allele:alleleset.elementSet()) {
				alleleStats[alleleCount][0] = allele;
				alleleStats[alleleCount++][1] = alleleset.count(allele);
			}
			
			//sort alleleStats by count, decreasing order
			if (nalleles > 1) {
				for (int i = 0; i < nalleles - 1; i++) for (int j = i + 1; j < nalleles; j++) {
					if (alleleStats[i][1] < alleleStats[j][1]) {
						int[] tmp = alleleStats[i];
						alleleStats[i] = alleleStats[j];
						alleleStats[j] = tmp;
					}
				}
			}
			
			alleleCounts.add(alleleStats);
		}
	}
	
	/**
	 * This function returns a copy of the component Haplotypes not references to them.
	 * @return	a deep copy of this cluster.
	 */
	public HaplotypeCluster copy() {
		return new HaplotypeCluster(new ArrayList<Haplotype>(hapList), score);
	}
	
	@Override
	public int compareTo(HaplotypeCluster cluster) {
		if (score > cluster.score) return -1;
		if (score < cluster.score) return 1;
		if (getSize() > cluster.getSize()) return -1;
		if (getSize() < cluster.getSize()) return 1;
		return 0;
	}
	
	/**
	 * The common haplotype of this cluster. The type is determined by the field ReturnHaplotype.
	 * @return a byte array representing the haplotype of this cluster
	 */
	public byte[] getHaplotype() {
		switch(ReturnHaplotype) {
		case unanimous:
			return getUnanimousHaplotype();
		case majority:
			return getMajorityHaplotype();
		case censored:
			return getCensoredMajorityHaplotype(0.05, 1);
		default:
			return getUnanimousHaplotype();
		}
	}
	
	/**
	 * The common haplotype of this cluster. The type is determined by the field ReturnHaplotype.
	 * @return	a String representing the haplotype of this cluster
	 */
	public String getHaplotypeAsString() {
		switch(ReturnHaplotype) {
		case unanimous:
			return byteHaplotypeAsString(getUnanimousHaplotype());
		case majority:
			return byteHaplotypeAsString(getMajorityHaplotype());
		default:
			return byteHaplotypeAsString(getUnanimousHaplotype());
		}
	}

	/**
	 * The majority haplotype of this cluster where each site is assigned value of the major allele within the cluster. 
	 * @return a byte array representing the haplotype of this cluster
	 */
	public byte[] getMajorityHaplotype() {
		int nsites = hapList.get(0).seqlen;
		byte[] hap = new byte[nsites];
		Arrays.fill(hap, N);
		
		for (int s = 0; s < nsites; s++) {
			//Get a count of all bytes not equal to N
			Multiset<Byte> byteset = HashMultiset.create(5);
			for (Haplotype h : hapList) {
				byte b = h.seq[s];
				if (b != N) byteset.add(h.seq[s]);
			}
			
			//if all bytes equal N, return N
			if (byteset.size() == 0) hap[s] = N;
			else {
				Iterator<Byte> it = byteset.elementSet().iterator();
				int maxcount = 0;
				Byte maxbyte = -1;
				boolean tie = false;
				while (it.hasNext()) {
					Byte thisbyte = it.next();
					int thiscount = byteset.count(thisbyte);
					if (thiscount > maxcount) {
						maxcount = thiscount;
						maxbyte = thisbyte;
						tie = false;
					} else if (thiscount == maxcount) {
						tie = true;
					}
				}
				if (tie) hap[s] = N;
				else hap[s] = maxbyte;
			}
		}
		
		return hap;
	}
	
	/**
	 * This returns the majority haplotype for sites with either minor allele frequency below maxMinorFreq or minor allele count less than or equal to maxMinorCount.
	 * Missing is returned at sites not meeting the stated criteria.
	 * @param maxMinorFreq	the maximum proportion of offtypes (alternate alleles) allowed before site is set to missing
	 * @param maxMinorCount	the maximum minor allele count allowed before site is set to missing	
	 * @return a byte array representing the haplotype of this cluster
	 */
	public byte[] getCensoredMajorityHaplotype(double maxMinorFreq, int maxMinorCount) {
		int nsites = hapList.get(0).seqlen;
		byte[] hap = new byte[nsites];
		Arrays.fill(hap, N);

		for (int s = 0; s < nsites; s++) {
			//Get a count of all bytes not equal to N
			Multiset<Byte> byteset = HashMultiset.create(5);
			for (Haplotype h : hapList) {
				byte b = h.seq[s];
				if (b != N) byteset.add(h.seq[s]);
			}

			//if all bytes are N, leave site as N
			if (byteset.size() > 0) {
				int majorcount = 0;
				int minorcount = 0;
				int totalcount = 0;
				Byte majorAllele = null;
				Byte minorAllele = null;
				Iterator<Byte> it = byteset.elementSet().iterator();
				
				while (it.hasNext()) {
					Byte thisbyte = it.next();
					int thiscount = byteset.count(thisbyte);
					totalcount += thiscount;
					if (thiscount > majorcount) {
						minorcount = majorcount;
						minorAllele = majorAllele;
						majorcount = thiscount;
						majorAllele = thisbyte;
					} else if (thiscount > minorcount) {
						minorcount = thiscount;
						minorAllele = thisbyte;
					}
				}
				double maf = ((double) minorcount) / totalcount;
				if (maf <= maxMinorFreq || minorcount <= maxMinorCount) {
					hap[s] = majorAllele.byteValue();
				}
			}
			
		}
		
		return hap;
	}
	
	/**
	 * The unanimous haplotype of this cluster. Only monomorphic sites are reported. Others will be unknown.
	 * @return a byte array representing the haplotype of this cluster
	 */
	public byte[] getUnanimousHaplotype() {
		int nsites = hapList.get(0).seqlen;
		byte[] hap = new byte[nsites];
		for (int s = 0; s < nsites; s++) {
			HashSet<Byte> byteset = new HashSet<Byte>();
			Multiset<Byte> multiByteset = HashMultiset.create();
			for (Haplotype haplo : hapList) {
				if (haplo.seq[s] != N) {
					byteset.add(haplo.seq[s]);
					multiByteset.add(haplo.seq[s]);
				}
			}
			
			int nNucleotidesObserved = byteset.size();
			if (nNucleotidesObserved == 1) { 
				hap[s] = Byte.valueOf(byteset.iterator().next());
			} else  {
				hap[s] = N;
			}
		}
		
		return hap;
	}
	
	/**
	 * Converts a byte sequence to a String
	 * @param hap an array of byte genotype or haplotype values
	 * @return the String representation
	 */
	public String byteHaplotypeAsString(byte[] hap) {
		StringBuilder sb = new StringBuilder();
		for (byte h:hap) sb.append(NucleotideAlignmentConstants.getNucleotideIUPAC(h));
		return sb.toString();
	}

	public int getNumberOfSites() {
		return hapList.get(0).seqlen;
	}
	
	public int[] listTaxaInCluster() {
		int ntaxa = hapList.size();
		int[] taxa = new int[ntaxa];
		int tcount = 0;
		for (Haplotype hap:hapList) taxa[tcount++] = hap.taxonIndex;
		Arrays.sort(taxa);
		return taxa;
	}
	
	public int countTaxaSharedWith(HaplotypeCluster hapCluster) {
		int[] otherTaxa = hapCluster.listTaxaInCluster();
		Arrays.sort(otherTaxa);
		int count = 0;
		
		for (Haplotype hap:hapList) {
			if (Arrays.binarySearch(otherTaxa, hap.taxonIndex) > -1) count++;
		}
		return count;
	}
	
	public boolean contains(Haplotype hap) {
		return hapList.contains(hap);
	}
	
	public int countHeterozygousSites() {
		countAllelesAtAllSites();
		int nhet = 0;
		int ntotal = getNumberOfSites();
		for (int s = 0; s < ntotal; s++) {
			int[][] alleles = alleleCounts.get(s);
			if (alleles.length > 0 && GenotypeTableUtils.isHeterozygous((byte) alleles[0][0])) nhet++;
			else if (alleles.length > 1 && alleles[1][1] > 1) nhet++;
		}
		return nhet;
	}
	
	@Override
	public String toString() {
		return getHaplotypeAsString();
	}
	
}
