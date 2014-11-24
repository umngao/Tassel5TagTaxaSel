package net.maizegenetics.dna.snp;

import java.util.Spliterator;
import java.util.function.Consumer;

import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.dna.map.Position;

public class ReferenceProbabilitySpliterator implements Spliterator<ReferenceProbabilitySpliterator.ReferenceProbabilityBySite> {
	private final GenotypeTable myGenotype;
	private final int numberOfTaxa;
	private int firstSite; 	//inclusive
	private int lastSite; 	//exclusive
	private final int minsize = 5;
	
	public ReferenceProbabilitySpliterator(GenotypeTable genotype, int firstSite, int lastSite) {
		if (firstSite < 0) throw new IllegalArgumentException("First site less than zero in ReferenceProbabilitySpliterator.");
		if (lastSite > genotype.numberOfSites()) throw new IllegalArgumentException("Last site too big in ReferenceProbabilitySpliterator.");
		myGenotype = genotype;
		this.firstSite = firstSite;
		this.lastSite = lastSite;
		numberOfTaxa = myGenotype.numberOfTaxa();
	}
	
	@Override
	public Spliterator<ReferenceProbabilitySpliterator.ReferenceProbabilityBySite> trySplit() {
		if (estimateSize() < minsize) return null;
		int mid = (lastSite - firstSite) /2;
		lastSite = mid;
		return new ReferenceProbabilitySpliterator(myGenotype, mid, lastSite);
	}

	@Override
	public long estimateSize() {
		return lastSite - firstSite;
	}

	@Override
	public int characteristics() {
		return Spliterator.IMMUTABLE | Spliterator.NONNULL | Spliterator.SIZED | Spliterator.SUBSIZED;
	}

	@Override
	public boolean tryAdvance(Consumer<? super ReferenceProbabilitySpliterator.ReferenceProbabilityBySite> action) {
		if (firstSite >= lastSite) return false;
		ReferenceProbabilitySpliterator.ReferenceProbabilityBySite info = new ReferenceProbabilitySpliterator.ReferenceProbabilityBySite();
		info.myPosition = myGenotype.positions().get(firstSite);
		float[] prob = new float[numberOfTaxa];
		for (int i = 0; i < numberOfTaxa; i++) prob[i] = myGenotype.referenceProbability(i, firstSite);
		info.myValues = AssociationUtils.convertFloatArrayToDouble(prob);
		action.accept(info);
		firstSite++;
		return true;
	}
	
	@Override
	public void forEachRemaining(Consumer<? super ReferenceProbabilitySpliterator.ReferenceProbabilityBySite> action) {
		for (int s = firstSite; s < lastSite; s++) {
			ReferenceProbabilitySpliterator.ReferenceProbabilityBySite info = new ReferenceProbabilitySpliterator.ReferenceProbabilityBySite();
			info.myPosition = myGenotype.positions().get(firstSite);
			float[] prob = new float[numberOfTaxa];
			for (int i = 0; i < numberOfTaxa; i++) prob[i] = myGenotype.referenceProbability(i, firstSite);
			info.myValues = AssociationUtils.convertFloatArrayToDouble(prob);
			action.accept(info);
		}
	}
	
	public class ReferenceProbabilityBySite {
		public Position myPosition;
		public double[] myValues;
		
		public ReferenceProbabilityBySite() {}
		public ReferenceProbabilityBySite(Position pos, double[] values) {
			myPosition = pos;
			myValues = values;
		}
	}
}


