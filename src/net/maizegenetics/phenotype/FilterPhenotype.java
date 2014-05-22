package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.TreeSet;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;


public class FilterPhenotype extends AbstractPhenotype {
	private final boolean myIsRowFilter;
	private final int[] myRowRedirect;
	private final int numberOfObservations;
	
	//private constructors
	
	private FilterPhenotype(AbstractPhenotype basePheno, String tablename, ArrayList<PhenotypeAttribute> retainedAttributes, ArrayList<ATTRIBUTE_TYPE> types, TaxaList taxaToRetain) {
		super(tablename, retainedAttributes, types);
		
		if (taxaToRetain == null) {
			myIsRowFilter = false;
			myRowRedirect = null;
			numberOfObservations = basePheno.getNumberOfObservations();
		} else if (taxaAttribute == null) {
			myLogger.error(String.format("Error in FilterPhenotype constructor: %s has no taxa", tablename));
			myIsRowFilter = false;
			myRowRedirect = null;
			numberOfObservations = basePheno.getNumberOfObservations();
		} else {
			
			myIsRowFilter = true;
			ArrayList<Taxon> retainedList = new ArrayList<>(taxaToRetain);
			Collections.sort(retainedList);
			int newObsNumber = 0;
			Taxon[] originalTaxa = (Taxon[]) taxaAttribute.getValues();
			int n = originalTaxa.length;
			int[] newObs = new int[n];
			for (int i = 0; i < n; i++) {
				if (Collections.binarySearch(retainedList, originalTaxa[i]) > -1) newObs[newObsNumber++] = i;
			}
			myRowRedirect = Arrays.copyOf(newObs, newObsNumber);
			numberOfObservations = newObsNumber;
		}
		
	}
	
	private FilterPhenotype(FilterPhenotype basePheno, String tablename, ArrayList<PhenotypeAttribute> retainedAttributes, ArrayList<ATTRIBUTE_TYPE> types, TaxaList taxaToRetain) {
		super(tablename, retainedAttributes, types);

		if (taxaToRetain == null) {
			myIsRowFilter = basePheno.myIsRowFilter;
			myRowRedirect = basePheno.myRowRedirect;
			numberOfObservations = basePheno.getNumberOfObservations();
		} else if (taxaAttribute == null) {
			myLogger.error(String.format("Error in FilterPhenotype constructor: %s has no taxa", tablename));
			myIsRowFilter = basePheno.myIsRowFilter;
			myRowRedirect = basePheno.myRowRedirect;
			numberOfObservations = basePheno.getNumberOfObservations();
		} else {
			myIsRowFilter = true;
			ArrayList<Taxon> retainedList = new ArrayList<>(taxaToRetain);
			Collections.sort(retainedList);
			int newObsNumber = 0;
			Taxon[] originalTaxa = (Taxon[]) taxaAttribute.getValues();
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
	
	public static Phenotype getInstance(Phenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, TaxaList taxaToRetain, String tablename) {
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>();
		for (PhenotypeAttribute attr : retainedAttributes) types.add(basePheno.getAttributeType(attr));
		if (tablename == null || tablename.trim().length() == 0) tablename = "filtered_" + basePheno.getName();

		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, tablename, retainedAttributes, types, taxaToRetain);
		else if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, tablename, retainedAttributes, types, taxaToRetain);
		else return null;
	}
	
	public static Phenotype getInstance(Phenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, String tablename) {
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>();
		for (PhenotypeAttribute attr : retainedAttributes) types.add(basePheno.getAttributeType(attr));
		if (tablename == null || tablename.trim().length() == 0) tablename = "filtered_" + basePheno.getName();
		
		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, tablename, retainedAttributes, types, null);
		else if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, tablename, retainedAttributes, types, null);
		else return null;
	}
	
	public static Phenotype getInstance(Phenotype basePheno, TaxaList retainedTaxa, String tablename) {
		if (tablename == null || tablename.trim().length() == 0) tablename = "filtered_" + basePheno.getName();
		
		if (basePheno instanceof CorePhenotype) {
			return new FilterPhenotype((CorePhenotype) basePheno, tablename, basePheno.getAttributeList(), basePheno.getAttributeTypeList(), retainedTaxa);
		}
		else if (basePheno instanceof FilterPhenotype) {
			return new FilterPhenotype((FilterPhenotype) basePheno, tablename, basePheno.getAttributeList(), basePheno.getAttributeTypeList(), retainedTaxa);
		}
		else return null;
	}
	
	public static Phenotype getInstanceRemoveTaxa(Phenotype basePheno, TaxaList taxaToRemove, String tablename) {
		ArrayList<Taxon> keepList = new ArrayList<Taxon>(basePheno.taxa());
		keepList.removeAll(taxaToRemove);
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		TaxaList keepTaxa = taxaBuilder.addAll(keepList).build();
		if (tablename == null || tablename.trim().length() == 0) tablename = "filtered_" + basePheno.getName();
		
		if (basePheno instanceof CorePhenotype) {
			return new FilterPhenotype((CorePhenotype) basePheno, tablename, basePheno.getAttributeList(), basePheno.getAttributeTypeList(), keepTaxa);
		}
		else if (basePheno instanceof FilterPhenotype) {
			return new FilterPhenotype((FilterPhenotype) basePheno, tablename, basePheno.getAttributeList(), basePheno.getAttributeTypeList(), keepTaxa);
		}
		else return null;
	}
	
	//public methods
	
	public int[] getRowRedirect() {
		return myRowRedirect;
	}
	
	@Override
	public int getNumberOfObservations() {
		return numberOfObservations;
	}

	@Override
	public Object getValue(int obs, int attrnum) {
		return super.getValue(myRowRedirect[obs], attrnum);
	}

	@Override
	public Object getValues(int attrnum) {
		return getAttribute(attrnum).getSubsetOfValues(myRowRedirect);
	}

	@Override
	public Object getValues(PhenotypeAttribute attr) {
		return attr.getSubsetOfValues(myRowRedirect);
	}

	@Override
	public TaxaList taxa() {
		TreeSet<Taxon> taxonset = new TreeSet<>();
		for (int ndx : myRowRedirect) {
			taxonset.add((Taxon) taxaAttribute.getValue(ndx));
		}
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		return taxaBuilder.addAll(taxonset).build();
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
	public Object[] getRow(int row) {
		return super.getRow(myRowRedirect[row]);
	}

	@Override
	public Object getValueAt(int row, int col) {
		return super.getValueAt(myRowRedirect[row], col);
	}

	@Override
	public int getRowCount() {
		return numberOfObservations;
	}


}
