package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Ints;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;

/**
 * @author pbradbury
 *
 */
public class CorePhenotype implements Phenotype {
	protected final List<PhenotypeAttribute> myAttributeList;
	protected final List<ATTRIBUTE_TYPE> myAttributeTypeList;
	protected final Multimap<ATTRIBUTE_TYPE, Integer> myAttributeTypeMap;
	protected final HashMap<String, Integer> myAttributeNameMap;
	protected final int numberOfAttributes;
	protected final int numberOfObservations;
	protected final String name;
	protected final TaxaAttribute myTaxaAttribute;
	
	CorePhenotype(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types, String name) {
		myAttributeList = new ArrayList<PhenotypeAttribute>(attributes);
		myAttributeTypeList = new ArrayList<ATTRIBUTE_TYPE>(types);
		this.name = name;
		myAttributeTypeMap = HashMultimap.create();
		int typeCount = 0;
		for (ATTRIBUTE_TYPE type:types) myAttributeTypeMap.put(type, typeCount++);
		numberOfAttributes = attributes.size();
		numberOfObservations = myAttributeList.get(0).size();
		
		Collection<Integer> taxaIndices = myAttributeTypeMap.get(ATTRIBUTE_TYPE.taxa);
		if (taxaIndices.size() == 1) {
			Integer Ndx = taxaIndices.iterator().next();
			myTaxaAttribute = (TaxaAttribute) myAttributeList.get(Ndx);
		} else myTaxaAttribute = null;
		
		myAttributeNameMap = new HashMap<>();
		int attrCount = 0;
		for (PhenotypeAttribute attr : attributes) myAttributeNameMap.put(attr.name(), attrCount++);
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
	public Object value(int obs, int attrnum) {
		return myAttributeList.get(attrnum).value(obs);
	}
	
	@Override
	public boolean isMissing(int obs, int attrnum) {
		return myAttributeList.get(attrnum).isMissing(obs);
	}

	@Override
	public PhenotypeAttribute attribute(int attrnum) {
		return myAttributeList.get(attrnum);
	}

	@Override
	public int indexOfAttribute(PhenotypeAttribute attribute) {
		return myAttributeList.indexOf(attribute);
	}

	@Override
	public List<PhenotypeAttribute> attributeListCopy() {
		return new ArrayList<PhenotypeAttribute>(myAttributeList);
	}

	@Override
	public List<ATTRIBUTE_TYPE> typeListCopy() {
		return new ArrayList<Phenotype.ATTRIBUTE_TYPE>(myAttributeTypeList);
	}

	@Override
	public TaxaList taxa() {
		int[] taxaIndices = attributeIndicesOfType(ATTRIBUTE_TYPE.taxa);
		if (taxaIndices.length == 1) {
			TaxaAttribute taxaAttr = (TaxaAttribute) myAttributeList.get(taxaIndices[0]);
			TreeSet<Taxon> taxaSet = new TreeSet<>();
			for (Taxon taxon : taxaAttr.allTaxa()) taxaSet.add(taxon);
			return new TaxaListBuilder().addAll(taxaSet).build();
		}
		return null;
	}

	@Override
	public int numberOfAttributes() {
		return numberOfAttributes;
	}

	@Override
	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type) {
		return attributeIndicesOfType(type).length;
	}

	@Override
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type) {
		return Ints.toArray(myAttributeTypeMap.get(type));
	}

	@Override
	public ATTRIBUTE_TYPE attributeType(int attrnum) {
		return myAttributeTypeList.get(attrnum);
	}

	@Override
	public int numberOfObservations() {
		return numberOfObservations;
	}

	@Override
	public String attributeName(int attrnum) {
		return myAttributeList.get(attrnum).name();
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public boolean hasTaxaAttribute() {
		return myTaxaAttribute != null;
	}

	@Override
	public TaxaAttribute taxaAttribute() {
		return myTaxaAttribute;
	}

	@Override
	public int attributeIndexForName(String name) {
		Integer Index = myAttributeNameMap.get(name);
		if (Index == null) return -1;
		return Index.intValue();
	}

	public static boolean areAttributeAndTypeListsCompatible(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types) {
		if (attributes.size() != types.size()) return false;
		boolean compatible = true;
		Iterator<ATTRIBUTE_TYPE> typeIter = types.iterator();
		for (PhenotypeAttribute attr : attributes) {
			ATTRIBUTE_TYPE myType = typeIter.next();
			if (!attr.isTypeCompatible(myType)) {
				compatible = false;
				break;
			}
		}
		return compatible;
	}
	
}
