package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.taxa.TaxaList;

public interface Phenotype {
	public enum ATTRIBUTE_TYPE {data, covariate, factor};

	public abstract Object getValue(int obs, int attrnum);

	public abstract PhenotypeAttribute getAttribute(int attrnum);

	/**
	 * @return	a shallow copy of the attribute list. 
	 */
	public abstract ArrayList<PhenotypeAttribute> getAttributeList();

	public abstract PhenotypeAttribute getAttributeOfType(ATTRIBUTE_TYPE type,
			int attrnum);

	public abstract List<PhenotypeAttribute> getAttributeListOfType(
			ATTRIBUTE_TYPE type);

	public abstract TaxaList taxa();

	public abstract int getNumberOfAttributes();

	public abstract int getNumberOfAttributesOfType(ATTRIBUTE_TYPE type);

	public abstract int getNumberOfObservations();

}