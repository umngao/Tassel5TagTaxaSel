package net.maizegenetics.phenotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
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
		if (source == SOURCE_TYPE.file) {
			try {
				BufferedReader br = new BufferedReader(new FileReader(filename));
				String topline = br.readLine();
				if (topline.toLowerCase().startsWith("<phenotype")) {
					String name = new File(filename).getName();
					if (name.endsWith(".txt")) name = name.substring(0, name.length() - 4);
					Phenotype myPhenotype = importPhenotypeFile(br, name);
					br.close();
					return myPhenotype;
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return null;
	}
	
	//private methods for importing phenotypes from different sources
	private Phenotype importPhenotypeFile(BufferedReader phenotypeReader, String name) {
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
				if (typeString[pheno].toLowerCase().startsWith("cov") || typeString[pheno].toLowerCase().equals("dat")) {
					float[] dataArray = new float[nObs];
					OpenBitSet missing = new OpenBitSet(nObs);
					int obsCount = 0;
					for (String[] inputLine : stringData) {
						try {
							dataArray[obsCount] = Float.parseFloat(inputLine[pheno]);
						} catch (NumberFormatException nfe) {
							dataArray[obsCount] = Float.NaN;
							missing.fastSet(pheno);
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
					attributes.add(new TaxaAttribute(taxa));
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
		
		return new CorePhenotype(attributes, types, name);
	}
	
}
