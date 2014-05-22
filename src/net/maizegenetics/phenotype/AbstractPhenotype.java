package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import com.google.common.collect.ArrayListMultimap;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;

public class AbstractPhenotype implements Phenotype, TableReport {
	protected final static Logger myLogger = Logger.getLogger(AbstractPhenotype.class);
	protected final ArrayListMultimap<ATTRIBUTE_TYPE, PhenotypeAttribute> attributeMultimap;
	protected final ArrayList<PhenotypeAttribute> myAttributeList;
	protected final ArrayList<ATTRIBUTE_TYPE> myAttributeTypes;
	protected final Map<PhenotypeAttribute, Integer> myAttributeIndex;
	protected final int numberOfAttributes;
	protected final int numberOfObservations;
	protected final String tableName;
	protected final PhenotypeAttribute taxaAttribute;

	protected AbstractPhenotype(String tableName, ArrayList<PhenotypeAttribute> attributes, ArrayList<ATTRIBUTE_TYPE> types) {
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
			taxaAttribute = attrs.get(0);
		} else {
			taxaAttribute = null;
		}
	}
	
	@Override
	public Object getValue(int obs, int attrnum) {
		return myAttributeList.get(attrnum).getValue(obs);
	}

	@Override
	public Object getValues(int attrnum) {
		return myAttributeList.get(attrnum).getValues();
	}

	@Override
	public Object getValues(PhenotypeAttribute attr) {
		return attr.getValues();
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
	public ArrayList<ATTRIBUTE_TYPE> getAttributeTypeList() {
		return new ArrayList<ATTRIBUTE_TYPE>(myAttributeTypes);
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
		if (taxaAttribute == null) return null;
		TreeSet<Taxon> taxaset = new TreeSet<>();
		Taxon[] taxa = (Taxon[]) taxaAttribute.getValues();
		for (Taxon taxon : taxa) taxaset.add(taxon);
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		return taxaBuilder.addAll(taxaset).build();
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
		ATTRIBUTE_TYPE oldtype = myAttributeTypes.get(attrnum);
		PhenotypeAttribute thisAttribute = myAttributeList.get(attrnum);
		myAttributeTypes.set(attrnum, type);
		attributeMultimap.remove(oldtype, thisAttribute);
		attributeMultimap.put(type, thisAttribute);
	}

	@Override
	public void setAttributeType(PhenotypeAttribute attribute, ATTRIBUTE_TYPE type) {
		int ndx = myAttributeIndex.get(attribute);
		ATTRIBUTE_TYPE oldtype = myAttributeTypes.get(ndx);
		myAttributeTypes.set(ndx, type);
		attributeMultimap.remove(oldtype, attribute);
		attributeMultimap.put(type, attribute);
	}

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
		return numberOfAttributes;
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
		Object[] rowValues = new Object[numberOfAttributes];
		for (int i = 0; i < numberOfAttributes; i++) rowValues[i] = getValueAt(row, i);
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
		return myAttributeList.get(col).getValue(row);
	}

	@Override
	public String getName() {
		return tableName;
	}

}