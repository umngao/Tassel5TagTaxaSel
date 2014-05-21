package net.maizegenetics.phenotype;

import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.ImmutableBiMap;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;

public class TaxaAttribute implements PhenotypeAttribute {
	TaxaList myTaxaList;
	private final int[] values;
	private final ImmutableBiMap<Taxon, Integer> taxaBimap;
	
	TaxaAttribute(List<Taxon> taxa) {
		TreeSet<Taxon> taxaSet = new TreeSet<>();
		taxaSet.addAll(taxa);
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		myTaxaList = taxaBuilder.addAll(taxaSet).build();
		ImmutableBiMap.Builder<Taxon, Integer> bimapBuilder = ImmutableBiMap.builder();
		int count = 0;
		for (Taxon taxon : taxaSet) bimapBuilder.put(taxon, count++);
		taxaBimap = bimapBuilder.build();
		
		int n = taxa.size();
		values = new int[n];
		for (int i = 0; i < n; i++) values[i] = taxaBimap.get(taxa.get(i));
	}
	
	@Override
	public Object getValue(int obs) {
		return taxaBimap.inverse().get(values[obs]);
	}

	@Override
	public Object getValues() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isMissing(int obs) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public BitSet getMissing() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getSize() {
		// TODO Auto-generated method stub
		return 0;
	}

}
