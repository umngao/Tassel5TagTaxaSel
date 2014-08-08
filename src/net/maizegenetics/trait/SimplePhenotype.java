package net.maizegenetics.trait;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

import java.util.List;

public class SimplePhenotype extends AbstractPhenotype {
	private DoubleMatrix2D data;
	
	public SimplePhenotype(TaxaList taxa, List<Trait> traits, DoubleMatrix2D data) {
		super(taxa, traits);
		this.data = data;
	}

	public SimplePhenotype(TaxaList taxa, List<Trait> traits, double[][] data) {
		super(taxa, traits);
		this.data = DoubleFactory2D.dense.make(data);
	}
	
	public SimplePhenotype(TaxaList taxa, List<Trait> traits) {
		super(taxa, traits);
		this.data = DoubleFactory2D.dense.make(getNumberOfTaxa(), getNumberOfTraits(), Double.NaN);
	}
	
	public double getData(int taxon, int trait) {
		return data.get(taxon, trait);
	}

	public double getData(Taxon taxon, Trait trait) {
		return getData(whichTaxon(taxon), whichTrait(trait));
	}

	public double[][] getData() {
		return data.toArray();
	}

	public void setData(int taxon, int trait, double value) {
		data.set(taxon, trait, value);
	}

	public void setData(Taxon taxon, Trait trait, double value) {
		data.set(whichTaxon(taxon), whichTrait(trait), value);
	}

}
