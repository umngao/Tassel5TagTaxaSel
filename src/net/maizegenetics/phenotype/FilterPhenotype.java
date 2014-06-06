package net.maizegenetics.phenotype;

import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Ints;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;

public class FilterPhenotype implements Phenotype, TableReport {
	private int[] myRowRedirect;
	private int numberOfObservations;
	private CorePhenotype basePhenotype;
	private String name;
	
	FilterPhenotype(CorePhenotype basePheno, List<Taxon> taxaToKeep, String name) {
		//TODO finish implementing taxaFilter
		this.name = name;
	}
	
	FilterPhenotype(FilterPhenotype basePheno, List<Taxon> taxaToKeep, String name) {
		//TODO finish implementing
		this.name = name;
	}
	
	FilterPhenotype(Phenotype basePheno, List<Taxon> taxaToKeep, String name) {
		//TODO finish implementing
		this.name = name;
	}
	
	//TableReport methods
	@Override
	public Object[] getTableColumnNames() {
		return basePhenotype.getTableColumnNames();
	}

	@Override
	public Object[][] getTableData() {
		Object[][] tableData = new Object[numberOfObservations][];
		for (int i = 0; i < numberOfObservations; i++) {
			tableData[i] = basePhenotype.getRow(myRowRedirect[i]);
		}
		return null;
	}

	@Override
	public String getTableTitle() {
		return name;
	}

	@Override
	public int getColumnCount() {
		return basePhenotype.getColumnCount();
	}

	@Override
	public int getRowCount() {
		return numberOfObservations;
	}

	@Override
	public int getElementCount() {
		return getRowCount() * getColumnCount();
	}

	@Override
	public Object[] getRow(int row) {
		return basePhenotype.getRow(myRowRedirect[row]);
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		int numberOfRows = end - start + 1;
		Object[][] tableData = new Object[numberOfRows][];
		for (int i = 0; i < numberOfRows; i++) {
			tableData[i] = basePhenotype.getRow(i + start);
		}
		return null;
	}

	@Override
	public Object getValueAt(int row, int col) {
		return basePhenotype.getValueAt(myRowRedirect[row], col);
	}

	//Phenotype methods
	
	@Override
	public Object value(int obs, int attrnum) {
		return basePhenotype.value(myRowRedirect[obs], attrnum);
	}

	@Override
	public boolean isMissing(int obs, int attrnum) {
		return basePhenotype.isMissing(myRowRedirect[obs], attrnum);
	}

	@Override
	public PhenotypeAttribute attribute(int attrnum) {
		return basePhenotype.attribute(attrnum).subset(myRowRedirect);
	}

	@Override
	public int indexOfAttribute(PhenotypeAttribute attribute) {
		return basePhenotype.indexOfAttribute(attribute);
	}

	@Override
	public List<PhenotypeAttribute> attributeListCopy() {
		return basePhenotype.attributeListCopy();
	}

	@Override
	public List<ATTRIBUTE_TYPE> typeListCopy() {
		return basePhenotype.typeListCopy();
	}

	@Override
	public int numberOfAttributes() {
		return basePhenotype.numberOfAttributes;
	}

	@Override
	public int numberOfObservations() {
		return numberOfObservations;
	}

	@Override
	public TaxaList taxa() {
		if (!basePhenotype.hasTaxaAttribute()) return null;
		Taxon[] taxaArray = basePhenotype.taxaAttribute().allTaxa();
		TreeSet<Taxon> taxaSet = new TreeSet<>();
		for (int ndx : myRowRedirect) taxaSet.add(taxaArray[ndx]);
		
		return new TaxaListBuilder().addAll(taxaSet).build();
	}

	@Override
	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type) {
		return basePhenotype.numberOfAttributesOfType(type);
	}

	@Override
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type) {
		return basePhenotype.attributeIndicesOfType(type);
	}

	@Override
	public ATTRIBUTE_TYPE attributeType(int attrnum) {
		return basePhenotype.attributeType(attrnum);
		
	}

	@Override
	public String attributeName(int attrnum) {
		return basePhenotype.attributeName(attrnum);
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public boolean hasTaxaAttribute() {
		return basePhenotype.hasTaxaAttribute();
	}

	@Override
	public TaxaAttribute taxaAttribute() {
		return (TaxaAttribute) basePhenotype.taxaAttribute().subset(myRowRedirect);
	}

}
