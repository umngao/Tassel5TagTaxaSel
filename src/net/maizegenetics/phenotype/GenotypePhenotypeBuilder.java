package net.maizegenetics.phenotype;

import net.maizegenetics.dna.snp.FilterGenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;

public class GenotypePhenotypeBuilder {
	Phenotype myPhenotype = null;
	GenotypeTable myGenotype = null;
	boolean isUnion = false;
	String myName = null;
	
	public GenotypePhenotypeBuilder phenotype(Phenotype thePhenotype) {
		myPhenotype = thePhenotype;
		return this;
	}
	
	public GenotypePhenotypeBuilder genotype(GenotypeTable theGenotype) {
		myGenotype = theGenotype;
		return this;
	}
	
	public GenotypePhenotypeBuilder union() {
		isUnion = true;
		return this;
	}
	
	public GenotypePhenotypeBuilder intersect() {
		isUnion = false;
		return this;
	}
	
	public GenotypePhenotypeBuilder name(String name) {
		myName = name;
		return this;
	}
	
	public GenotypePhenotype build() {
		if (myPhenotype == null) throw new IllegalArgumentException("Error: no phenotype data set was specified.");
		if (myGenotype == null) throw new IllegalArgumentException("Error: no genotype data set was specified.");
		if (myName == null) {
			myName = myPhenotype.name() + "_with_genotypes";
		}
		
		if (!isUnion) {
			TaxaList commonTaxa = TaxaListUtils.getCommonTaxa(myGenotype.taxa(), myPhenotype.taxa());
			if (myGenotype.taxa().numberOfTaxa() > commonTaxa.numberOfTaxa()) {
				myGenotype = FilterGenotypeTableBuilder.getInstance(myGenotype).taxaToKeep(commonTaxa).build();
			}
			if (myPhenotype.taxa().numberOfTaxa() > commonTaxa.numberOfTaxa()) {
				myPhenotype = new PhenotypeBuilder().filterPhenotype(myPhenotype).keepTaxa(commonTaxa).build();
			}
		}
		
		return new GenotypePhenotype(myGenotype, myPhenotype, myName);
	}
}
