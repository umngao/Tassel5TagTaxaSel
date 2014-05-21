package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;

public class FilterPhenotype implements Phenotype, TableReport {
	private final boolean myIsRowFilter;
	private final int[] myRowRedirect;
	private final ArrayList<PhenotypeAttribute> myAttributes; 
	private final ArrayList<ATTRIBUTE_TYPE> myAttributeTypes;
	private final ArrayListMultimap<ATTRIBUTE_TYPE, PhenotypeAttribute> attributeMultimap;
	private final TaxaAttribute myTaxa;
	private final int numberOfAttributes;
	private final int numberOfObservations;
	
	private final String myTableTitle = "Filtered phenotype";
	
	//private constructors
	
	private FilterPhenotype(CorePhenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, TaxaList taxaToRetain) {
		if (retainedAttributes == null) {
			myAttributes = basePheno.getAttributeList();
		} else {
			myAttributes = new ArrayList<PhenotypeAttribute>(retainedAttributes);
		}
		
		numberOfAttributes = myAttributes.size();
		myAttributeTypes = new ArrayList<>();
		attributeMultimap = ArrayListMultimap.create();
		for (PhenotypeAttribute attr : myAttributes) {
			ATTRIBUTE_TYPE myType = basePheno.getAttributeType(attr);
			myAttributeTypes.add(myType);
			attributeMultimap.put(myType, attr);
		}
		
		List<PhenotypeAttribute> taxaAttrList = attributeMultimap.get(ATTRIBUTE_TYPE.taxa);
		if (taxaAttrList.size() != 1) myTaxa = null;
		else myTaxa = (TaxaAttribute) taxaAttrList.get(0);

		if (taxaToRetain == null) {
			myIsRowFilter = false;
			myRowRedirect = null;
			numberOfObservations = basePheno.getNumberOfObservations();
		} else {
			myIsRowFilter = true;
			ArrayList<Taxon> retainedList = new ArrayList<>(taxaToRetain);
			Collections.sort(retainedList);
			int newObsNumber = 0;
			Taxon[] originalTaxa = (Taxon[]) myTaxa.getValues();
			int n = originalTaxa.length;
			int[] newObs = new int[n];
			for (int i = 0; i < n; i++) {
				if (Collections.binarySearch(retainedList, originalTaxa[i]) > -1) newObs[newObsNumber++] = i;
			}
			myRowRedirect = Arrays.copyOf(newObs, newObsNumber);
			numberOfObservations = newObsNumber;
		}
		
		
	}
	
	private FilterPhenotype(FilterPhenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, TaxaList taxaToRetain) {
		//TODO optimize for FilterPhenotype
		if (retainedAttributes == null) {
			myAttributes = basePheno.getAttributeList();
		} else {
			myAttributes = new ArrayList<PhenotypeAttribute>(retainedAttributes);
		}
		
		numberOfAttributes = myAttributes.size();
		myAttributeTypes = new ArrayList<>();
		attributeMultimap = ArrayListMultimap.create();
		for (PhenotypeAttribute attr : myAttributes) {
			ATTRIBUTE_TYPE myType = basePheno.getAttributeType(attr);
			myAttributeTypes.add(myType);
			attributeMultimap.put(myType, attr);
		}
		
		List<PhenotypeAttribute> taxaAttrList = attributeMultimap.get(ATTRIBUTE_TYPE.taxa);
		if (taxaAttrList.size() != 1) myTaxa = null;
		else myTaxa = (TaxaAttribute) taxaAttrList.get(0);

		if (taxaToRetain == null) {
			myIsRowFilter = false;
			myRowRedirect = null;
			numberOfObservations = basePheno.getNumberOfObservations();
		} else {
			myIsRowFilter = true;
			ArrayList<Taxon> retainedList = new ArrayList<>(taxaToRetain);
			Collections.sort(retainedList);
			int newObsNumber = 0;
			Taxon[] originalTaxa = (Taxon[]) myTaxa.getValues();
			int n = originalTaxa.length;
			int[] newObs = new int[n];
			for (int i = 0; i < n; i++) {
				if (Collections.binarySearch(retainedList, originalTaxa[i]) > -1) newObs[newObsNumber++] = i;
			}
			myRowRedirect = Arrays.copyOf(newObs, newObsNumber);
			numberOfObservations = newObsNumber;
		}
		
	}
	
	//public static getInstance methods
	
	public static Phenotype getInstance(Phenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, TaxaList taxaToRetain) {
		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, retainedAttributes, taxaToRetain);
		else if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, retainedAttributes, taxaToRetain);
		else return null;
	}
	
	public static Phenotype getInstance(Phenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes) {
		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, retainedAttributes, null);
		else if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, retainedAttributes, null);
		else return null;
	}
	
	public static Phenotype getInstance(Phenotype basePheno, TaxaList retainedTaxa) {
		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, null, retainedTaxa);
		else if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, null, retainedTaxa);
		else return null;
	}
	
	public static Phenotype getInstanceRemoveTaxa(Phenotype basePheno, TaxaList taxaToRemove) {
		ArrayList<Taxon> keepList = new ArrayList<Taxon>(basePheno.taxa());
		keepList.removeAll(taxaToRemove);
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		TaxaList keepTaxa = taxaBuilder.addAll(keepList).build();
		
		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, null, keepTaxa);
		else if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, null, keepTaxa);
		else return null;
	}
	

	@Override
	public int getAttributeIndex(PhenotypeAttribute attr) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public ATTRIBUTE_TYPE getAttributeType(int attrnum) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ATTRIBUTE_TYPE getAttributeType(PhenotypeAttribute attribute) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setAttributeType(int attrnum, ATTRIBUTE_TYPE type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setAttributeType(PhenotypeAttribute attribute,
			ATTRIBUTE_TYPE type) {
		// TODO Auto-generated method stub
		
	}

	//methods to override ------------------------------------
	@Override
	public Object getValue(int obs, int attrnum) {
		if (myIsRowFilter) {
			return myAttributes.get(attrnum).getValue(myRowRedirect[obs]);
		} else {
			return myAttributes.get(attrnum).getValue(obs);
		}
	}

	@Override
	public PhenotypeAttribute getAttribute(int attrnum) {
		if (myIsRowFilter) {
			return myAttributes.get(attrnum).getSubset(myRowRedirect);	
		} else return myAttributes.get(attrnum);
	}

	@Override
	public ArrayList<PhenotypeAttribute> getAttributeList() {
		return new ArrayList<PhenotypeAttribute>(myAttributes);
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
		if (myTaxa == null) return null;
		
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		Taxon[] taxa;
		if (myIsRowFilter) {
			taxa = (Taxon[]) myTaxa.getSubsetOfValues(myRowRedirect);
		} else {
			taxa = (Taxon[]) myTaxa.getValues();
		}
		
		return taxaBuilder.addAll(taxa).build();
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
	public Object[] getTableColumnNames() {
		String[] names = new String[numberOfAttributes];
		for (int a = 0; a < numberOfAttributes; a++) {
			names[a] = myAttributes.get(a).getName();
		}
		return names;
	}

	@Override
	public Object[][] getTableData() {
		int nrows = numberOfObservations;
		int ncols = numberOfAttributes;
		Object[][] resultTable = new Object[nrows][ncols];
		for (int r = 0; r < nrows; r++) {
			for (int c = 0; c < ncols; c++) {
				resultTable[r][c] = getValueAt(r,c);
			}
		}
		return resultTable;
	}

	@Override
	public String getTableTitle() {
		return myTableTitle;
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
	public Object[] getRow(int row) {
		Object[] values = new Object[numberOfAttributes];
		int obs;
		if (myIsRowFilter) obs = myRowRedirect[row];
		else obs = row;
		for (int i = 0; i < numberOfAttributes; i++) {
			values[i] = myAttributes.get(i).getValue(obs);
		}
		return values;
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		int nrows = end - start + 1;
		Object[][] resultTable = new Object[nrows][];
		for (int r = start; r <= end; r++) resultTable[r] = getRow(r);
		return resultTable;
	}

	@Override
	public Object getValueAt(int row, int col) {
		int obs;
		if (myIsRowFilter) obs = myRowRedirect[row];
		else obs = row;
		return myAttributes.get(col).getValue(obs);
	}

	@Override
	public int getElementCount() {
		return getRowCount() * getColumnCount();
	}

}
