package net.maizegenetics.phenotype;

import java.util.List;


import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class TaxaAttribute implements PhenotypeAttribute {
	private final List<Taxon> taxaList;
	private final int numberOfTaxa;
	
	TaxaAttribute(List<Taxon> taxa) {
		taxaList = taxa;
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
