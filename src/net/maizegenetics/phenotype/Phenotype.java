package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.taxa.TaxaList;

/**
 * Phenotype represents phenotype data as a two dimensional matrix. Rows are observations. Columns are attributes.
 * 
 * @author pbradbury
 *
 */
public interface Phenotype {
	public enum ATTRIBUTE_TYPE {data, covariate, factor, taxa};

	/**
	 * @param obs	an observation number
	 * @param attrnum	the index or column number of the attribute
	 * @return	the attribute value for the observation
	 */
	public abstract Object getValue(int obs, int attrnum);
	
	/**
	 * @param attrnum	the index or column number of the attribute
	 * @return	a 1 dimensional array of values for this attribute in observation order
	 * The type of array returned is a function of the attribute class. The array will be a primitive type where feasible.
	 * A CategoricalAttribute returns an array of int. A NumericAttribute returns an array of float. A TaxaAttribute returns an array of Taxon.
	 */
	public abstract Object getValues(int attrnum);
	
	/**
	 * @param attr	an attribute
	 * @return	a 1 dimensional array of values for this attribute in observation order
	 * The type of array returned is a function of the attribute class. The array will be a primitive type where feasible.
	 * A CategoricalAttribute returns an array of int. A NumericAttribute returns an array of float. A TaxaAttribute returns an array of Taxon.
	 */
	public abstract Object getValues(PhenotypeAttribute attr);

	/**
	 * @param attrnum	the index or column number of the attribute
	 * @return	the attribute for column attrnum
	 */
	public abstract PhenotypeAttribute getAttribute(int attrnum);
	
	/**
	 * @param attr	an attribute
	 * @return	the index or column number for attribute attr
	 */
	public abstract int getAttributeIndex(PhenotypeAttribute attr);

	/**
	 * @return	a shallow copy of the attribute list 
	 */
	public abstract ArrayList<PhenotypeAttribute> getAttributeList();
	
	/**
	 * @return a shallow copy of the attribute type list
	 */
	public abstract ArrayList<ATTRIBUTE_TYPE> getAttributeTypeList();

	public abstract PhenotypeAttribute getAttributeOfType(ATTRIBUTE_TYPE type,
			int attrnum);

	public abstract List<PhenotypeAttribute> getAttributeListOfType(
			ATTRIBUTE_TYPE type);

	public abstract TaxaList taxa();

	public abstract int getNumberOfAttributes();

	public abstract int getNumberOfAttributesOfType(ATTRIBUTE_TYPE type);
	
	public abstract ATTRIBUTE_TYPE getAttributeType(int attrnum);
	
	public abstract ATTRIBUTE_TYPE getAttributeType(PhenotypeAttribute attribute);
	
	public abstract void setAttributeType(int attrnum, ATTRIBUTE_TYPE type);
	
	public abstract void setAttributeType(PhenotypeAttribute attribute, ATTRIBUTE_TYPE type);

	public abstract int getNumberOfObservations();
	
	public abstract String getName();

}