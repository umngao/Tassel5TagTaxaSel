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
	
	public boolean isMissing(int obs, int attrnum);
	
	public PhenotypeAttribute attribute(int attrnum);

	public int numberOfAttributes();
	
	public int numberOfObservations();

	public TaxaList taxa();

	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type);
	
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type);
	
	public ATTRIBUTE_TYPE attributeType(int attrnum);
	
	public void setAttributeType(int attrnum, ATTRIBUTE_TYPE type);
	
	public String attributeName(int attrnum);
	
	public String name();
	
	public boolean hasTaxaAttribute();
	
	public TaxaAttribute taxaAttribute();

}