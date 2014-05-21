package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.TableReport;

import com.google.common.collect.ArrayListMultimap;

public class CorePhenotype implements TableReport, Phenotype {
	private ArrayListMultimap<ATTRIBUTE_TYPE, PhenotypeAttribute> attributeMultimap;
	private ArrayList<PhenotypeAttribute> myAttributeList;
	private TaxaList myTaxaList;
	private int numberOfAttributes;
	private int numberOfObservations;
	private String tableName;
	
	CorePhenotype(String tableName, ArrayList<PhenotypeAttribute> attributes, ArrayList<ATTRIBUTE_TYPE> types) {
		myAttributeList = attributes;
		numberOfAttributes = myAttributeList.size();
		numberOfObservations = attributes.get(0).getSize();
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
