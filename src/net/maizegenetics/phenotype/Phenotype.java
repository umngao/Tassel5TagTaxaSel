package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;

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
	public Object value(int obs, int attrnum);
	
	/**
	 * @param obs	an observation number
	 * @param attrnum	the index or column number of the attribute
	 * @return	true if the phenotype value of attrnum is missing for obs
	 */
	public boolean isMissing(int obs, int attrnum);
	
	/**
	 * @param attrnum	the index or column number of the attribute
	 * @return	the PhenotypeAttribute 
	 */
	public PhenotypeAttribute attribute(int attrnum);

	/**
	 * @param attribute	a PhenotypeAttribute
	 * @return	the index of attribute in this Phenotype, -1 if the attribute is not contained in this Phenotype
	 */
	public int indexOfAttribute(PhenotypeAttribute attribute);
	
	/**
	 * @return	the number of attributes or columns in the Phenotype
	 */
	public int numberOfAttributes();
	
	/**
	 * @return	the number of observations or rows in this phenotype
	 */
	public int numberOfObservations();

	/**
	 * @return	the unique set of taxa in this Phenotype. 
	 * Each taxon will be in the TaxaList only once though there may be more than one observation of the taxon in the Phenotype.
	 */
	public TaxaList taxa();

	/**
	 * @param type	an attribute type
	 * @return	the number of attributes of this type
	 */
	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type);
	
	/**
	 * @param type	an attribute type
	 * @return	an array of the indices of all the attributes of this type
	 */
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type);
	
	/**
	 * @param attrnum	the index or column number of an attribute
	 * @return	the type of the attribute
	 */
	public ATTRIBUTE_TYPE attributeType(int attrnum);
	
	/**
	 * @param attrnum	the index or column number of the attribute
	 * @return	the name of the attribute
	 */
	public String attributeName(int attrnum);
	
	/**
	 * @return	the name of this Phenotype
	 */
	public String name();
	
	/**
	 * @return	true, if this Phenotype type has one and only one TaxaAttribute
	 */
	public boolean hasTaxaAttribute();
	
	/**
	 * @return	this Phenotype's TaxaAttribute, which is the Taxa represented by the observations in this Phenotype
	 */
	public TaxaAttribute taxaAttribute();
	
	/**
	 * @return	a shallow copy of the attribute list for this Phenotype
	 */
	public List<PhenotypeAttribute> attributeListCopy();
	
	/**
	 * @return	a shallow copy of the attribute type list for this Phenotype
	 */
	public List<ATTRIBUTE_TYPE> typeListCopy();

}