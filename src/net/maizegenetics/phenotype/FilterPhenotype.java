package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;

public class FilterPhenotype implements Phenotype, TableReport {
	private CorePhenotype basePhenotype;
	private boolean myIsTaxaFilter = false;
	private boolean myIsAttributeFilter = false;
	private int[] myTaxaRedirect = null;
	private int[] myAttributeRedirect = null;
	Map<ATTRIBUTE_TYPE, int[]> myAttributeTypeRedirect = null;
	private String myTableTitle = "Filtered phenotype";
	
	private FilterPhenotype(CorePhenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, TaxaList taxaToRetain) {
		basePhenotype = basePheno;
		
		if (taxaToRetain != null) {
			myIsTaxaFilter = true;
			TaxaList commonTaxa = TaxaListUtils.getCommonTaxa(basePhenotype.taxa(), taxaToRetain);
			int ntaxa = commonTaxa.size();
			myTaxaRedirect = new int[ntaxa];
			int tcount = 0;
			for (Taxon taxon : commonTaxa) {
				myTaxaRedirect[tcount++] = basePhenotype.taxa().indexOf(taxon);
			}
		}
		
		if (retainedAttributes != null) {
			myIsAttributeFilter = true;
			int nAttr = retainedAttributes.size();
			myAttributeRedirect = new int[nAttr];
			int acount = 0;
			for (PhenotypeAttribute attr : retainedAttributes) {
				myAttributeRedirect[acount++] = basePhenotype.getAttributeList().indexOf(attr);
			}
		}
	}
	
	public static Phenotype getInstance(CorePhenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes, TaxaList taxaToRetain) {
		return new FilterPhenotype(basePheno, retainedAttributes, taxaToRetain);
	}
	
	public static Phenotype getInstance(CorePhenotype basePheno, ArrayList<PhenotypeAttribute> retainedAttributes) {
		return new FilterPhenotype(basePheno, retainedAttributes, null);
	}
	
	public static Phenotype getInstance(CorePhenotype basePheno, TaxaList retainedTaxa) {
		return new FilterPhenotype(basePheno, null, retainedTaxa);
	}
	
	public static Phenotype getInstanceRemoveTaxa(CorePhenotype basePheno, TaxaList taxaToRemove) {
		ArrayList<Taxon> keepList = new ArrayList<Taxon>(basePheno.taxa());
		keepList.removeAll(taxaToRemove);
		TaxaListBuilder taxaBuilder = new TaxaListBuilder();
		taxaBuilder.addAll(keepList);
		return new FilterPhenotype(basePheno, null, taxaBuilder.build());
	}
	

	//methods to override ------------------------------------
	@Override
	public Object getValue(int obs, int attrnum) {
		if (myIsTaxaFilter && myIsAttributeFilter) {
			return basePhenotype.getValue(myTaxaRedirect[obs], myAttributeRedirect[attrnum]);
		} else if (myIsAttributeFilter) {
			return basePhenotype.getValue(obs, myAttributeRedirect[attrnum]);
		} else {
			return basePhenotype.getValue(myTaxaRedirect[obs], attrnum);
		}
	}

	@Override
	public PhenotypeAttribute getAttribute(int attrnum) {
		if (myIsAttributeFilter) {
			return basePhenotype.getAttribute(myAttributeRedirect[attrnum]);	
		} else return basePhenotype.getAttribute(attrnum);
	}

	@Override
	public ArrayList<PhenotypeAttribute> getAttributeList() {
		if (myIsAttributeFilter) {
			ArrayList<PhenotypeAttribute> attrList = new ArrayList<PhenotypeAttribute>();
			for (int attrnum : myAttributeRedirect) attrList.add(basePhenotype.getAttribute(myAttributeRedirect[attrnum]));
			return attrList;	
		} else return basePhenotype.getAttributeList();
	}

	@Override
	public PhenotypeAttribute getAttributeOfType(ATTRIBUTE_TYPE type, int attrnum) {
		if (myIsAttributeFilter) {
			return basePhenotype.getAttribute(myAttributeTypeRedirect.get(type)[attrnum]);
		} else {
			return basePhenotype.getAttributeOfType(type, attrnum);
		}
	}

	@Override
	public List<PhenotypeAttribute> getAttributeListOfType(ATTRIBUTE_TYPE type) {
		if (myIsAttributeFilter) {
			ArrayList<PhenotypeAttribute> attrList = new ArrayList<PhenotypeAttribute>();
			for (int attrnum : myAttributeTypeRedirect.get(type)) attrList.add(basePhenotype.getAttribute(attrnum));
			return attrList;
		} else {
			return basePhenotype.getAttributeListOfType(type);
		}

	}

	@Override
	public TaxaList taxa() {
		if (myIsTaxaFilter) {
			TaxaListBuilder taxaBuilder = new TaxaListBuilder();
			for (int t : myTaxaRedirect) taxaBuilder.add(basePhenotype.taxa().get(t));
			return taxaBuilder.build();
		} else {
			return basePhenotype.taxa();
		}
	}

	@Override
	public int getNumberOfAttributes() {
		if (myIsAttributeFilter) {
			return myAttributeRedirect.length;
		} else {
			return basePhenotype.getNumberOfAttributes();
		}
	}

	@Override
	public int getNumberOfAttributesOfType(ATTRIBUTE_TYPE type) {
		if (myIsAttributeFilter) {
			return myAttributeTypeRedirect.get(type).length;
		} else return basePhenotype.getNumberOfAttributesOfType(type);
	}

	@Override
	public int getNumberOfObservations() {
		if (myIsTaxaFilter) {
			return myTaxaRedirect.length;
		} else return basePhenotype.getNumberOfObservations();
	}

	@Override
	public Object[] getTableColumnNames() {
		if (myIsAttributeFilter) {
			int nAttr = myAttributeRedirect.length;
			Object[] names = new Object[nAttr];
			for (int a = 0; a < nAttr; a++) {
				names[a] = basePhenotype.getAttribute(a).getName();
			}
			return names;
		} else return basePhenotype.getTableColumnNames();
	}

	@Override
	public Object[][] getTableData() {
		int nrows = myTaxaRedirect.length;
		int ncols = myAttributeRedirect.length + 1;
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
		if (myIsAttributeFilter) return myAttributeRedirect.length + 1;
		else return basePhenotype.getColumnCount();
	}

	@Override
	public int getRowCount() {
		 if (myIsTaxaFilter) return myTaxaRedirect.length;
		 else return basePhenotype.getRowCount();
	}

	@Override
	public Object[] getRow(int row) {
		int ncols = getColumnCount();
		Object[] rowData = new Object[ncols];
		for (int i = 0; i < ncols; i++) rowData[i] = getValueAt(row, i);
		return basePhenotype.getRow(row);
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		int nrows = end - start + 1;
		int ncols = myAttributeRedirect.length + 1;
		Object[][] resultTable = new Object[nrows][ncols];
		for (int r = 0; r < nrows ; r++) {
			for (int c = 0; c < ncols; c++) {
				resultTable[r][c] = getValueAt(r + start,c);
			}
		}
		return resultTable;
	}

	@Override
	public Object getValueAt(int row, int col) {
		if (myIsTaxaFilter && myIsAttributeFilter) {
			if (col == 0) return basePhenotype.taxa().get(myTaxaRedirect[row]);
			return basePhenotype.getValue(myTaxaRedirect[row], myAttributeRedirect[col - 1]);
		} else if (myIsAttributeFilter) {
			if (col == 0) return basePhenotype.taxa().get(row);
			return basePhenotype.getValue(row, myAttributeRedirect[col - 1]);
		} else {
			if (col == 0) return basePhenotype.taxa().get(myTaxaRedirect[row]);
			return basePhenotype.getValue(myTaxaRedirect[row], col - 1);
		}
	}

	@Override
	public int getElementCount() {
		return getRowCount() * getColumnCount();
	}

}
