package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 * @author pbradbury
 *
 */
public class TaxaAttribute implements PhenotypeAttribute {
	private final ArrayList<Taxon> taxaList;
	private final int numberOfTaxa;
	
	TaxaAttribute(List<Taxon> taxa) {
		if (taxa instanceof ArrayList) taxaList = (ArrayList<Taxon>) taxa;
		else taxaList = new ArrayList<>(taxa);
		numberOfTaxa = taxa.size();
	}
	
	@Override
	public Object getValue(int obs) {
		return taxaList.get(obs);
	}

	@Override
	public Object getValues() {
		Taxon[] taxaArray = new Taxon[numberOfTaxa];
		return taxaList.toArray(taxaArray);
	}

	@Override
	public Object getSubsetOfValues(int[] obs) {
		int n = obs.length;
		Taxon[] taxaArray = new Taxon[n];
		for (int i = 0; i < n;  i++) taxaArray[i] = taxaList.get(obs[i]);
		return taxaArray;
	}

	@Override
	public PhenotypeAttribute getSubset(int[] obs) {
		int n = obs.length;
		List<Taxon> subset = new ArrayList<Taxon>();
		for (int i = 0; i < n; i++) subset.add(taxaList.get(obs[i]));
		return null;
	}

	@Override
	public boolean isMissing(int obs) {
		return false;
	}

	@Override
	public BitSet getMissing() {
		return new OpenBitSet(numberOfTaxa);
	}

	@Override
	public String getName() {
		return "Taxa";
	}

	@Override
	public int getSize() {
		return numberOfTaxa;
	}

}
