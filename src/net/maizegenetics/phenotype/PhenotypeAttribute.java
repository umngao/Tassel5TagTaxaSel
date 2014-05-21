package net.maizegenetics.phenotype;

import net.maizegenetics.util.BitSet;

public interface PhenotypeAttribute {

	/**
	 * @param obs	the observation number
	 * @return	the value for this observation
	 */
	Object getValue(int obs);
	
	/**
	 * The return value will typically be a primitive array whose type depends on the sub class
	 * @return	the values of this Attribute for all observations, in order by observation number
	 */
	Object getValues();
	
	/**
	 * @param obs	an array of observation numbers 
	 * @return	the values of this attribute for the observations in obs in the same order
	 */
	Object getSubsetOfValues(int[] obs);
	
	/**
	 * @param obs	an array of observation numbers
	 * @return	an attribute equivalent to this one but with a subset of observations specified by obs
	 */
	PhenotypeAttribute getSubset(int[] obs);
	
	/**
	 * @param obs	the observation number
	 * @return	if the value of the observation is missing, true, otherwise false.
	 */
	boolean isMissing(int obs);
	
	/**
	 * @return an array whose elements are true for each missing observation, false for observations with valid values.
	 */
	BitSet getMissing();
	
	/**
	 * @return	the name of this Attribute
	 */
	String getName();
	
	/**
	 * @return	the number of observations in this Attribute
	 */
	int getSize();
}
