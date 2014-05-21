package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;

import com.google.common.collect.ArrayListMultimap;

public class CorePhenotype implements TableReport, Phenotype {
	private final ArrayListMultimap<ATTRIBUTE_TYPE, PhenotypeAttribute> attributeMultimap;
	private final ArrayList<PhenotypeAttribute> myAttributeList;
	private final ArrayList<ATTRIBUTE_TYPE> myAttributeTypes;
	private final Map<PhenotypeAttribute, Integer> myAttributeIndex;
	private final TaxaList myTaxaList;
	private final int numberOfAttributes;
	private final int numberOfObservations;
	private final String tableName;
	
	CorePhenotype(String tableName, ArrayList<PhenotypeAttribute> attributes, ArrayList<ATTRIBUTE_TYPE> types) {
		this.tableName = tableName;
		myAttributeList = attributes;
		myAttributeTypes = types;
		numberOfAttributes = myAttributeList.size();
		numberOfObservations = attributes.get(0).getSize();
		
		attributeMultimap = ArrayListMultimap.create();
		myAttributeIndex = new HashMap<PhenotypeAttribute, Integer>();
		
		for (int i = 0; i < numberOfAttributes; i++) {
			attributeMultimap.put(types.get(i), attributes.get(i));
			myAttributeIndex.put(attributes.get(i), i);
		}
		
		List<PhenotypeAttribute> attrs = attributeMultimap.get(ATTRIBUTE_TYPE.taxa);
		if (attrs.size() == 1) {
			TaxaAttribute taxaAttr = (TaxaAttribute) attrs.get(0);
			TreeSet<Taxon> taxaset = new TreeSet<>();
			Taxon[] taxa = (Taxon[]) taxaAttr.getValues();
			for (Taxon taxon : taxa) taxaset.add(taxon);
			TaxaListBuilder taxaBuilder = new TaxaListBuilder();
			myTaxaList = taxaBuilder.addAll(taxaset).build();
			
		} else myTaxaList = null;
	}
	
	@Override
	public Object getValue(int obs, int attrnum) {
		return myAttributeList.get(attrnum).getValue(obs);
	}
	
	@Override
	public PhenotypeAttribute getAttribute(int attrnum) {
		return myAttributeList.get(attrnum);
	}
	
	@Override
	public int getAttributeIndex(PhenotypeAttribute attr) {
		return myAttributeIndex.get(attr);
	}

	@Override
	public ArrayList<PhenotypeAttribute> getAttributeList() {
		return new ArrayList<PhenotypeAttribute>(myAttributeList);
	}
	
	@Override
	public PhenotypeAttribute getAttributeOfType(ATTRIBUTE_TYPE type, int attrnum) {
		return attributeMultimap.get(type).get(attrnum);
	}
	
	@Override
	public List<PhenotypeAttribute> getAttributeListOfType(ATTRIBUTE_TYPE type) {
		return new ArrayList<PhenotypeAttribute>(attributeMultimap.get(type));
	}
	
	@Override
	public TaxaList taxa() {
		return myTaxaList;
	}
	
	@Override
	public int getNumberOfAttributes() {
		return numberOfAttributes;
	}
	
	@Override
	public int getNumberOfAttributesOfType(ATTRIBUTE_TYPE type) {
		return attributeMultimap.get(type).size();
	}
	
	@Override
	public int getNumberOfObservations() {
		return numberOfObservations;
	}

	@Override
	public ATTRIBUTE_TYPE getAttributeType(int attrnum) {
		return myAttributeTypes.get(attrnum);
	}

	@Override
	public ATTRIBUTE_TYPE getAttributeType(PhenotypeAttribute attribute) {
		return myAttributeTypes.get(getAttributeIndex(attribute));
	}

	@Override
	public void setAttributeType(int attrnum, ATTRIBUTE_TYPE type) {
		myAttributeTypes.set(attrnum, type);
	}

	@Override
	public void setAttributeType(PhenotypeAttribute attribute,
			ATTRIBUTE_TYPE type) {
		myAttributeTypes.set(getAttributeIndex(attribute), type);
	}

	// TableReport functions
	@Override
	public Object[] getTableColumnNames() {
		String[] columnNames = new String[numberOfAttributes];
		int count = 0;
		for (PhenotypeAttribute attr : myAttributeList) {
			columnNames[count++] = attr.getName();
		}
		return columnNames;
	}
	
	@Override
	public Object[][] getTableData() {
		return getTableData(0, numberOfObservations - 1);
	}
	
	@Override
	public String getTableTitle() {
		return tableName;
	}
	
	@Override
	public int getColumnCount() {
		//the first column is the taxa names
		return numberOfAttributes + 1;
	}
	
	@Override
	public int getRowCount() {
		return numberOfObservations;
	}
	
	@Override
	public int getElementCount() {
		return getColumnCount() * getRowCount();
	}
	
	@Override
	public Object[] getRow(int row) {
		int n = numberOfAttributes + 1;
		Object[] rowValues = new Object[numberOfAttributes + 1];
		for (int i = 0; i < n; i++) rowValues[i] = getValueAt(row, i);
		return rowValues;
	}
	
	@Override
	public Object[][] getTableData(int start, int end) {
		int n = end - start + 1;
		Object[][] data = new Object[n][];
		for (int i = 0; i < n; i++) data[i] = getRow(start + i);
		return data;
	}
	
	@Override
	public Object getValueAt(int row, int col) {
		if (col == 0) return myTaxaList.get(row);
		return myAttributeList.get(col - 1).getValue(row);
	}
	
}
