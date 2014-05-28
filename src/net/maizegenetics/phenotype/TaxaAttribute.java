package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
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
	private final String name;
	
	private static final List<ATTRIBUTE_TYPE> myAllowedTypes;
	static{
		myAllowedTypes = new ArrayList<ATTRIBUTE_TYPE>();
		myAllowedTypes.add(ATTRIBUTE_TYPE.taxa);
	}

	TaxaAttribute(List<Taxon> taxa, String name) {
		this.name = name;
		if (taxa instanceof ArrayList) taxaList = (ArrayList<Taxon>) taxa;
		else taxaList = new ArrayList<>(taxa);
		numberOfTaxa = taxa.size();
	}
	
	TaxaAttribute(List<Taxon> taxa) {
		this(taxa, "Taxa");
	}
	
	public Taxon[] allTaxa() {
		return taxaList.toArray(new Taxon[numberOfTaxa]);
	}
	
	public Taxon taxon(int obs) {
		return taxaList.get(obs);
	}
	
	@Override
	public Object value(int obs) {
		return taxaList.get(obs);
	}

	@Override
	public Object allValues() {
		Taxon[] taxaArray = new Taxon[numberOfTaxa];
		return taxaList.toArray(taxaArray);
	}

	@Override
	public PhenotypeAttribute subset(int[] obs) {
		int n = obs.length;
		ArrayList<Taxon> subset = new ArrayList<Taxon>();
		for (int i = 0; i < n; i++) subset.add(taxaList.get(obs[i]));
		return new TaxaAttribute(subset);
	}

	@Override
	public boolean isMissing(int obs) {
		return false;
	}

	@Override
	public BitSet missing() {
		return new OpenBitSet(numberOfTaxa);
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public int size() {
		return numberOfTaxa;
	}

	@Override
	public List<ATTRIBUTE_TYPE> getCompatibleTypes() {
		return myAllowedTypes;
	}

	@Override
	public boolean isTypeCompatible(ATTRIBUTE_TYPE type) {
		return myAllowedTypes.contains(type);
	}

}
