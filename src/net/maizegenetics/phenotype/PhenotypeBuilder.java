package net.maizegenetics.phenotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.OpenBitSet;

/**
 * @author Peter Bradbury
 *
 */
public class PhenotypeBuilder {
	private Logger myLogger = Logger.getLogger(PhenotypeBuilder.class);
	private enum SOURCE_TYPE{file, phenotype, list, join};
	private SOURCE_TYPE source;
	private String filename;
	private Phenotype basePhenotype;
	private List<Phenotype> phenotypesToJoin;
	private List<Taxon> taxaToKeep = null;
	private List<Taxon> taxaToRemove = null;
	private List<PhenotypeAttribute> attributeList = null;
	private List<ATTRIBUTE_TYPE> attributeTypeList = null;
  	private int[] indexOfAttributesToKeep = null;
	private HashMap<ATTRIBUTE_TYPE, Integer> attributeChangeMap = new HashMap<ATTRIBUTE_TYPE, Integer>();
	private boolean isUnionJoin;
	private boolean isFilterable = false;
	private String phenotypeName = "Phenotype";
	public PhenotypeBuilder() {
		
	}
	
	/**
	 * @param basePhenotype	a base Phenotype to be filtered
	 * @return	a builder that can be used to filter a Phenotype
	 * Filtering means to generate a new Phenotype containing only a subset of the attributes and/or taxa of the original.
	 */
	public static PhenotypeBuilder getFilterableInstance(Phenotype basePhenotype) {
		return new PhenotypeBuilder(basePhenotype);
	}
	
	private PhenotypeBuilder(Phenotype basePhenotype) {
		this.basePhenotype = basePhenotype;
		source = SOURCE_TYPE.phenotype;
		isFilterable = true;
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
	 * @param phenotypes	a List of Phenotypes to be joined
	 * @return	a Phenotype builder that will build a join of the list of Phenotypes
	 * A union join returns a Phenotype containing any taxon present in at least one of the Phenotypes to be joined.
	 * An intersect join returns a Phenotype containing only taxa present in all of the Phenotypes to be joined.
	 * The type of join method should be specified using the intersect() or union() method. If no join method is specified, an intersect join will be performed.
	 */
	public PhenotypeBuilder joinPhenotypes(List<Phenotype> phenotypes) {
		phenotypesToJoin = phenotypes;
		isUnionJoin = false;
		source = SOURCE_TYPE.join;
		return this;
	}
	
	/**
	 * @return	a builder that will perform an intersect join if given a list of Phenotypes to join
	 */
	public PhenotypeBuilder intersect() {
		isUnionJoin = false;
		return this;
	}
	
	/**
	 * @return	a builder that will perform a union join if given a list of Phenotypes to join
	 */
	public PhenotypeBuilder union() {
		isUnionJoin = true;
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
	
	public PhenotypeBuilder assignName(String name) {
		phenotypeName = name;
		return this;
	}
	
	/**
	 * @param taxaToKeep	a list of taxa to be kept from the base Phenotype
	 * @return	a PhenotypeBuilder that will return a FilterPhenotype with taxa in the taxaToKeep list
	 * Only taxa that are in both taxaToKeep and the base Phenotype will be included in the Phenotype that is built.
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder keepTaxa(List<Taxon> taxaToKeep) {
		if (!isFilterable) notFilterable();
		this.taxaToKeep = taxaToKeep; 
		return this;
	}
	
	/**
	 * @param taxaToRemove	a list of taxa to removed from the base Phenotype
	 * @return	a PhenotypeBuilder that will return a Phenotype with taxa from the supplied list excluded
	 * Any taxa in taxaToRemove but not in the base Phenotype will be ignored.
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder removeTaxa(List<Taxon> taxaToRemove)  {
		if (!isFilterable) notFilterable();
		this.taxaToRemove = taxaToRemove;
		return this;
	}
	
	/**
	 * @param attributesToKeep	a list of the attributes to be kept
	 * @return	a PhenotypeBuilder that will return a new Phenotype with only the attributes in the supplied list
	 * Only attributes in both attributesToKeep and the base Phenotype will be included in the Phenotype that is built.
	 * The order of the attributes in the new Phenotype will match that in attributesToKeep.
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder keepAttributes(List<PhenotypeAttribute> attributesToKeep) {
		if (!isFilterable) notFilterable();
		attributeList = attributesToKeep;
		return this;
	}
	
	/**
	 * @param indexOfAttributes	the column numbers of the attributes in the base Phenotype to be included in the newly built Phenotype
	 * @return	a PhenotypeBuilder that will build a Phenotype with the specified attributes
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder keepAttributes(int[] indexOfAttributes) {
		if (!isFilterable) notFilterable();
		this.indexOfAttributesToKeep = indexOfAttributes;
		return this;
	}
	
	/**
	 * @param attributeIndex	the numeric index (column number) of an attribute in the base Phenotype
	 * @param type	the new type for that attribute
	 * @return	a PhenotypeBuilder that will build a phenotype with the changed attribute type
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder changeAttributeType(int attributeIndex, ATTRIBUTE_TYPE type) {
		if (!isFilterable) notFilterable();
		attributeChangeMap.put(type, attributeIndex );
		return this;
	}
	
	/**
	 * @param attributeTypes	a list of attribute types for the attributes to be built
	 * @return	a PhenotypeBuilder that will build a Phenotype that will have this list of types
	 * The order of types must be the same as the order of attributes as supplied by the keepAttributes methods if used or in the base Phenotype if the attribute list is not changed.
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder typesOfRetainedAttributes(List<ATTRIBUTE_TYPE> attributeTypes) {
		if (!isFilterable) notFilterable();
		attributeTypeList = attributeTypes;
		return this;
	}
	
	/**
	 * @return a new Phenotype built with the supplied parameters
	 */
	public Phenotype build() {
		if (source == SOURCE_TYPE.file) {
			try {
				BufferedReader br = new BufferedReader(new FileReader(filename));
				String topline = br.readLine();
				if (phenotypeName.equals("Phenotype")) {
					phenotypeName = new File(filename).getName();
					if (phenotypeName.endsWith(".txt")) phenotypeName = phenotypeName.substring(0, phenotypeName.length() - 4);
				}
				Phenotype myPhenotype;
				if (topline.toLowerCase().startsWith("<phenotype")) {
					myPhenotype = importPhenotypeFile(br);
				} else {
					myPhenotype = importTraitFile(br, topline);
				}
				br.close();
				return myPhenotype;
				
			} catch (IOException e) {
				e.printStackTrace();
				return null;
			}
		} else if (source == SOURCE_TYPE.phenotype) {
			return filterBasePhenotype();
		} else if (source == SOURCE_TYPE.list) {
			return createPhenotypeFromLists();
		} else if (source == SOURCE_TYPE.join) {
			return joinPhenotypes();
		}
		return null;
	}
	
	//private methods  ------------------------------------------------------
	private void notFilterable() {
		throw new java.lang.IllegalStateException("Phenotype Builder error: applied a filter method to a non-filterable instance.");
	}
	
	private Phenotype importPhenotypeFile(BufferedReader phenotypeReader) {
		Pattern whiteSpace = Pattern.compile("\\s+");
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>();
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>();
		
		//assumes the first line has been read to determine that this is indeed a Phenotype file
		try {
			String[] typeString = whiteSpace.split(phenotypeReader.readLine());
			String[] phenoNames = whiteSpace.split(phenotypeReader.readLine());
			int nPheno = typeString.length;
			ArrayList<String[]> stringData = new ArrayList<String[]>();
			String inputStr;
			while ((inputStr = phenotypeReader.readLine()) != null) {
				stringData.add(whiteSpace.split(inputStr));
			}
			
			int nObs = stringData.size();
			for (int pheno = 0; pheno < nPheno; pheno++) {
				if (typeString[pheno].toLowerCase().startsWith("cov") || typeString[pheno].toLowerCase().startsWith("dat")) {
					float[] dataArray = new float[nObs];
					OpenBitSet missing = new OpenBitSet(nObs);
					int obsCount = 0;
					for (String[] inputLine : stringData) {
						try {
							dataArray[obsCount] = Float.parseFloat(inputLine[pheno]);
						} catch (NumberFormatException nfe) {
							dataArray[obsCount] = Float.NaN;
							missing.fastSet(obsCount);
						}
						obsCount++;
					}
					attributes.add(new NumericAttribute(phenoNames[pheno], dataArray, missing));
					if (typeString[pheno].toLowerCase().startsWith("cov")) types.add(ATTRIBUTE_TYPE.covariate);
					else types.add(ATTRIBUTE_TYPE.data);
				} else if (typeString[pheno].toLowerCase().startsWith("tax")) {
					ArrayList<Taxon> taxa = new ArrayList<>();
					for (String[] inputLine : stringData) {
						taxa.add(new Taxon(inputLine[pheno]));
					}
					attributes.add(new TaxaAttribute(taxa, phenoNames[pheno]));
					types.add(ATTRIBUTE_TYPE.taxa);
				} else if (typeString[pheno].toLowerCase().startsWith("fac")) {
					String[] labelArray = new String[nObs];
					int obsCount = 0;
					for (String[] inputLine : stringData) {
						labelArray[obsCount++] = inputLine[pheno];
					}
					attributes.add(new CategoricalAttribute(phenoNames[pheno], labelArray));
					types.add(ATTRIBUTE_TYPE.factor);
				}
			}
			
		} catch (IOException e) {
			myLogger.error("Error reading phenotype file.");
			e.printStackTrace();
		}
		
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private Phenotype importTraitFile(BufferedReader phenotypeReader, String firstLine) {
		//TODO implement
		return null;
	}
	
	private Phenotype filterBasePhenotype() {
		//TODO implement
		return null;
	}
	
	private Phenotype createPhenotypeFromLists() {
		if (attributeList.size() != attributeTypeList.size()) throw new IllegalArgumentException("Error building Phenotype: attribute list size not equal to type list size.");
		Iterator<ATTRIBUTE_TYPE> typeIter = attributeTypeList.iterator();
		for (PhenotypeAttribute attr : attributeList) {
			ATTRIBUTE_TYPE type = typeIter.next();
			if (!attr.isTypeCompatible(type)) {
				throw new IllegalArgumentException("Error building Phenotype: types not compatible with attributes.");
			}
		}
		return new CorePhenotype(attributeList, attributeTypeList, phenotypeName);
	}
	
	private Phenotype joinPhenotypes() {
		//TODO implement
		return null;
	}
}
