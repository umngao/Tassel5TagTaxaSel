package net.maizegenetics.phenotype;

import java.util.List;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.TableReport;

public class CorePhenotype implements Phenotype, TableReport {
	protected final List<PhenotypeAttribute> myAttributeList;
	protected final List<ATTRIBUTE_TYPE> myAttributeTypeList;
	protected final int numberOfAttributes;
	protected final int numberOfObservations;
	protected final String name;
	
	public CorePhenotype() {
		myAttributeList = null;
		myAttributeTypeList = null;
		numberOfAttributes = 0;
		numberOfObservations = 0;
		name = "Phenotype";
	}
	
	//TableReport methods
	
	@Override
	public Object[] getTableColumnNames() {
		String[] names = new String[numberOfAttributes];
		int ptr = 0;
		for (PhenotypeAttribute attr : myAttributeList) names[ptr++] = attr.name();
		return names;
	}

	@Override
	public Object[][] getTableData() {
		Object[][] tableData = new Object[numberOfObservations][numberOfAttributes];
		int ptr = 0;
		for (PhenotypeAttribute attr : myAttributeList) {
			for (int i = 0; i < numberOfObservations; i++) tableData[i][ptr] = attr.value(i);
			ptr++;
		}
		return tableData;
	}

	@Override
	public String getTableTitle() {
		return name;
	}

	@Override
	public int getColumnCount() {
		return numberOfAttributes;
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
		Object[] rowData = new Object[numberOfAttributes];
		int ptr = 0;
		for (PhenotypeAttribute attr : myAttributeList) rowData[ptr++] = attr.value(row);
		return rowData;
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		int numberOfRows = end - start + 1;
		Object[][] tableData = new Object[numberOfRows][numberOfAttributes];
		int ptr = 0;
		for (PhenotypeAttribute attr : myAttributeList) {
			for (int i = 0; i < numberOfRows; i++) tableData[i][ptr] = attr.value(i + start);
			ptr++;
		}
		return tableData;
	}

	@Override
	public Object getValueAt(int row, int col) {
		return myAttributeList.get(col).value(row);
	}

	//Phenotype methods
	
	@Override
	public Object getValue(int obs, int attrnum) {
		return myAttributeList.get(attrnum).value(obs);
	}

	@Override
	public boolean isMissing(int obs, int attrnum) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public TaxaList taxa() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int numberOfAttributes() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ATTRIBUTE_TYPE attributeType(int attrnum) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void attributeType(int attrnum, ATTRIBUTE_TYPE type) {
		// TODO Auto-generated method stub

	}

	@Override
	public int numberOfObservations() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String attributeName(int attrnum) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String name() {
		// TODO Auto-generated method stub
		return null;
	}

}
