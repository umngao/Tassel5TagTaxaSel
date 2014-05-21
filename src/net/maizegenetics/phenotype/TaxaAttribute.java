package net.maizegenetics.phenotype;

import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.ImmutableBiMap;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class TaxaAttribute implements PhenotypeAttribute {
	private final TaxaList myTaxaList;
	private final int[] values;
	private final ImmutableBiMap<Taxon, Integer> taxaBimap;
	private final int numberOfValues;
	
	TaxaAttribute(List<Taxon> taxa) {
		TreeSet<Taxon> taxaSet = new TreeSet<>();
		taxaSet.addAll(taxa);
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		myTaxaList = taxaBuilder.addAll(taxaSet).build();
		ImmutableBiMap.Builder<Taxon, Integer> bimapBuilder = ImmutableBiMap.builder();
		int count = 0;
		for (Taxon taxon : taxaSet) bimapBuilder.put(taxon, count++);
		taxaBimap = bimapBuilder.build();
		
		numberOfValues = taxa.size();
		values = new int[numberOfValues];
		for (int i = 0; i < numberOfValues; i++) values[i] = taxaBimap.get(taxa.get(i));
		
	}
	
	@Override
	public Object getValue(int obs) {
		return taxaBimap.inverse().get(values[obs]);
	}

	@Override
	public Object getValues() {
		Taxon[] taxa = new Taxon[numberOfValues];
		for (int i = 0; i < numberOfValues; i++) taxa[i] = taxaBimap.inverse().get(i);
		return taxa;
	}

	@Override
	public boolean isMissing(int obs) {
		return false;
	}

	@Override
	public BitSet getMissing() {
		return new OpenBitSet(numberOfValues);
	}

	@Override
	public String getName() {
		return "Taxa";
	}

	@Override
	public int getSize() {
		return numberOfValues;
	}

}
