package net.maizegenetics.phenotype;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.Taxon;

/**
 * @author Peter Bradbury
 *
 */
public class PhenotypeBuilder {
	private enum SOURCE_TYPE{file, phenotype, list};
	private SOURCE_TYPE source;
	private String filename;
	private Phenotype basePhenotype;
	private List<Taxon> taxaToKeep = null;
	private List<Taxon> taxaToRemove = null;
	private List<PhenotypeAttribute> attributeList = null;
	private List<ATTRIBUTE_TYPE> attributeTypeList = null;
  	private int[] indexOfAttributesToKeep = null;
	private HashMap<ATTRIBUTE_TYPE, Integer> attributeChangeMap = new HashMap<ATTRIBUTE_TYPE, Integer>();
	
	public PhenotypeBuilder() {
		
	}
	 
	/**
	 * @param filename	the name of a file containing phenotype data to be imported
	 * @return	a PhenotypeBuilder that will import a file
	 */
	public PhenotypeBuilder fromFile(String filename) {
		this.filename = filename;
		source = SOURCE_TYPE.file;
		return this;
	}
	
	/**
	 * @param basePhenotype	the Phenotype to be modified or just copied
	 * @return	this PhenotypeBuilder that will be used to filter the base Phenotype
	 */
	public PhenotypeBuilder fromPhenotype(Phenotype basePhenotype) {
		this.basePhenotype = basePhenotype;
		source = SOURCE_TYPE.phenotype;
		return this;
	}
	
	/**
	 * @param attributes	a list of attributes
	 * @param types	a list of types matching the attribute list
	 * @return	a PhenotypeBuilder that will build using these lists
	 * The attribute and type lists must be the same size and the types must be compatible with the attributes
	 */
	public PhenotypeBuilder fromList(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types) {
		attributeList = attributes;
		attributeTypeList = types;
		source = SOURCE_TYPE.list;
		return this;
	}
	
	/**
	 * @param taxaToKeep	a list of taxa to be kept from the base Phenotype
	 * @return	a PhenotypeBuilder that will return a FilterPhenotype with taxa in the taxaToKeep list
	 * Only taxa that are in both taxaToKeep and the base Phenotype will be included in the Phenotype that is built.
	 */
	public PhenotypeBuilder keepTaxa(List<Taxon> taxaToKeep) {
		this.taxaToKeep = taxaToKeep; 
		return this;
	}
	
	/**
	 * @param taxaToRemove	a list of taxa to removed from the base Phenotype
	 * @return	a PhenotypeBuilder that will return a Phenotype with taxa from the supplied list excluded
	 * Any taxa in taxaToRemove but not in the base Phenotype will be ignored.
	 */
	public PhenotypeBuilder removeTaxa(List<Taxon> taxaToRemove)  {
		this.taxaToRemove = taxaToRemove;
		return this;
	}
	
	/**
	 * @param attributesToKeep	a list of the attributes to be kept
	 * @return	a PhenotypeBuilder that will return a new Phenotype with only the attributes in the supplied list
	 * Only attributes in both attributesToKeep and the base Phenotype will be included in the Phenotype that is built.
	 * The order of the attributes in the new Phenotype will match that in attributesToKeep.
	 */
	public PhenotypeBuilder keepAttributes(List<PhenotypeAttribute> attributesToKeep) {
		attributeList = attributesToKeep;
		return this;
	}
	
	/**
	 * @param indexOfAttributes	the column numbers of the attributes in the base Phenotype to be included in the newly built Phenotype
	 * @return	a PhenotypeBuilder that will build a Phenotype with the specified attributes
	 */
	public PhenotypeBuilder keepAttributes(int[] indexOfAttributes) {
		this.indexOfAttributesToKeep = indexOfAttributes;
		return this;
	}
	
	/**
	 * @param attributeIndex	the numeric index (column number) of an attribute in the base Phenotype
	 * @param type	the new type for that attribute
	 * @return	a PhenotypeBuilder that will build a phenotype with the changed attribute type
	 */
	public PhenotypeBuilder changeAttributeType(int attributeIndex, ATTRIBUTE_TYPE type) {
		attributeChangeMap.put(type, attributeIndex );
		return this;
	}
	
	/**
	 * @param attributeTypes	a list of attribute types for the attributes to be built
	 * @return	a PhenotypeBuilder that will build a Phenotype that will have this list of types
	 * The order of types must be the same as the order of attributes as supplied by the keepAttributes methods if used or in the base Phenotype if the attribute list is not changed.
	 */
	public PhenotypeBuilder typesOfRetainedAttributes(List<ATTRIBUTE_TYPE> attributeTypes) {
		attributeTypeList = attributeTypes;
		return this;
	}
	
	/**
	 * @return a new Phenotype built with the supplied parameters
	 */
	public Phenotype build() {
		//TODO implement
		return null;
	}
	
}
