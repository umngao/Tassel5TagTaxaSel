package net.maizegenetics.phenotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
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
	private boolean isUnionJoin = false;
	private boolean isFilterable = false;
	private boolean isConcatenate = false;
	private String phenotypeName = "Phenotype";
	private List<PhenotypeAttribute> attributeList = null;
	private List<ATTRIBUTE_TYPE> attributeTypeList = null;
	
	//filter criteria
	private List<Taxon> taxaToKeep = null;
	private List<Taxon> taxaToRemove = null;
  	private int[] indexOfAttributesToKeep = null;
	private HashMap<PhenotypeAttribute, ATTRIBUTE_TYPE> attributeChangeMap = new HashMap<PhenotypeAttribute, ATTRIBUTE_TYPE>();
	
	public PhenotypeBuilder() {
		
	}
	
	/**
	 * @param basePhenotype	the base Phenotype that will be filtered or copied
	 * @return	a filterable PhenotypeBuilder
	 * This is the only method that returns a filterable PhenotypeBuilder
	 */
	public PhenotypeBuilder filterPhenotype(Phenotype basePhenotype) {
		this.basePhenotype = basePhenotype;
		source = SOURCE_TYPE.phenotype;
		phenotypeName = "filtered_" + basePhenotype.name();
		isFilterable = true;
		return this;
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
		source = SOURCE_TYPE.join;
		StringBuilder sb = new StringBuilder();
		for (Phenotype pheno : phenotypes) {
			if (sb.length() > 0) sb.append(" + ");
			sb.append(pheno.name());
		}
		phenotypeName = sb.toString();
		return this;
	}
	
	/**
	 * @param pheno1	a Phenotype
	 * @param pheno2	a second Phenotype to be merged with pheno1
	 * @return	a Phenotype builder that will merge pheno1 and pheno2
	 */
	public PhenotypeBuilder joinPhenotypes(Phenotype pheno1, Phenotype pheno2) {
		List<Phenotype> phenoList = new ArrayList<Phenotype>();
		phenoList.add(pheno1);
		phenoList.add(pheno2);
		return joinPhenotypes(phenoList);
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
	 * @return	a builder that will perform concatenate the phenotypes to be joined
	 */
	public PhenotypeBuilder concatenate() {
		isConcatenate = true;
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
		taxaToRemove = null;
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
		taxaToKeep = null;
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
		indexOfAttributesToKeep = null;
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
		attributeList = null;
		return this;
	}
	
	/**
	 * @param attributeIndex	the numeric index (column number) of an attribute in the base Phenotype
	 * @param type	the new type for that attribute
	 * @return	a PhenotypeBuilder that will build a phenotype with the changed attribute type
	 * This function can only be applied to a filterable instance.
	 */
	public PhenotypeBuilder changeAttributeType(PhenotypeAttribute attribute, ATTRIBUTE_TYPE type) {
		if (!isFilterable) notFilterable();
		attributeChangeMap.put(attribute, type);
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
				File phenotypeFile = new File(filename);
				BufferedReader br = new BufferedReader(new FileReader(phenotypeFile));
				String topline = br.readLine();
				if (phenotypeName.equals("Phenotype")) {
					phenotypeName = new File(filename).getName();
					if (phenotypeName.endsWith(".txt")) phenotypeName = phenotypeName.substring(0, phenotypeName.length() - 4);
				}
				Phenotype myPhenotype;
				if (topline.toLowerCase().startsWith("<phenotype")) {
					myPhenotype = importPhenotypeFile(phenotypeFile);
				} else {
					myPhenotype = importTraits(phenotypeFile);
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
		throw new java.lang.IllegalStateException("Phenotype Builder error: applied a filter method to a non-filterable builder.");
	}
	
	private Phenotype importPhenotypeFile(File phenotypeFile) throws IOException {
		Pattern whiteSpace = Pattern.compile("\\s+");
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>();
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>();

		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));
		phenotypeReader.readLine();  //assumes the first line has been read to determine that this is indeed a Phenotype file
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
		phenotypeReader.close();

		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private Phenotype importTraits(File phenotypeFile) throws IOException {
		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));

		int numberOfDataLines = 0;
		String inputline = phenotypeReader.readLine();
		ArrayList<String> headerLines = new ArrayList<>();
		boolean isFactor = false;
		boolean isCovariate = false;
		boolean hasHeaders = false;
		boolean isTrait = false;
		String[] traitnames = new String[0];
		while (inputline != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				numberOfDataLines++;
			} else if (inputline.toLowerCase().startsWith("<trai")) {
				isTrait = true;
				String[] splitLine = inputline.split("[<>\\s]+");
				traitnames = Arrays.copyOfRange(splitLine, 2, splitLine.length);
			} else if (inputline.toLowerCase().startsWith("<cov")) {
				isCovariate = true;
			} else if (inputline.toLowerCase().startsWith("<fac")) {
				isFactor = true;
			} else if (inputline.toLowerCase().startsWith("<head")) {
				hasHeaders = true;
				headerLines.add(inputline);
			}
			
			inputline = phenotypeReader.readLine();
		}
		phenotypeReader.close();
		
		if (hasHeaders) {
			processTraitsAndFactors(phenotypeFile, traitnames, numberOfDataLines, isCovariate, headerLines);
		} else if (isFactor) {
			return processFactors(phenotypeFile, traitnames, numberOfDataLines);
		} else if (isTrait) {
			return processTraits(phenotypeFile, traitnames, numberOfDataLines, isCovariate);
		}
		throw new IllegalArgumentException("Unrecognized format for a phenotype.");
	}
	
	private Phenotype processTraits(File phenotypeFile, String[] traitnames, int numberOfDataLines, boolean isCovariate) throws IOException {
		int ntraits = traitnames.length;
		int nattributes = ntraits + 1;
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(nattributes);
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(nattributes);
		ArrayList<float[]> traitValues = new ArrayList<>(ntraits);
		ArrayList<BitSet> missingList = new ArrayList<>(ntraits);
		ArrayList<Taxon> taxaList = new ArrayList<>();
		
		types.add(ATTRIBUTE_TYPE.taxa);
		ATTRIBUTE_TYPE myAttributeType;
		if (isCovariate) myAttributeType = ATTRIBUTE_TYPE.covariate;
		else myAttributeType = ATTRIBUTE_TYPE.data;
		for (int i = 0; i < ntraits; i++) {
			traitValues.add(new float[numberOfDataLines]);
			missingList.add(new OpenBitSet(numberOfDataLines));
			types.add(myAttributeType);
		}
		
		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));
		int dataCount = 0;
		int lineCount = 1;
		String inputline;
		while ((inputline = phenotypeReader.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				String[] values = inputline.split("\\s+");
				if (values.length != ntraits + 1) {
					String msg = String.format("Incorrect number of values in line %d of %s", lineCount, phenotypeFile.getName());
					phenotypeReader.close();
					throw new IllegalArgumentException(msg);
				}
				taxaList.add(new Taxon(values[0]));
				for (int i = 0; i < ntraits; i++) {
					float val;
					try {
						val = Float.parseFloat(values[i + 1]);
					} catch (NumberFormatException e) {
						val = Float.NaN;
						missingList.get(i).fastSet(dataCount);
					}
					traitValues.get(i)[dataCount] = val;
				}
				dataCount++;
			}
			
			lineCount++;
		}
		
		phenotypeReader.close();
		
		attributes.add(new TaxaAttribute(taxaList));
		for (int i = 0; i < ntraits; i++) attributes.add(new NumericAttribute(traitnames[i], traitValues.get(i), missingList.get(i)));
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private Phenotype processFactors(File phenotypeFile, String[] traitnames, int numberOfDataLines) throws IOException {
		int ntraits = traitnames.length;
		int nattributes = ntraits + 1;
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(nattributes);
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(nattributes);
		ArrayList<String[]> traitValues = new ArrayList<>(ntraits);
		ArrayList<Taxon> taxaList = new ArrayList<>();
		
		types.add(ATTRIBUTE_TYPE.taxa);
		ATTRIBUTE_TYPE myAttributeType = ATTRIBUTE_TYPE.factor;
		for (int i = 0; i < ntraits; i++) {
			traitValues.add(new String[numberOfDataLines]);
			types.add(myAttributeType);
		}
		
		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));
		int dataCount = 0;
		int lineCount = 1;
		String inputline;
		while ((inputline = phenotypeReader.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				String[] values = inputline.split("\\s+");
				if (values.length != ntraits + 1) {
					String msg = String.format("Incorrect number of values in line %d of %s", lineCount, phenotypeFile.getName());
					phenotypeReader.close();
					throw new IllegalArgumentException(msg);
				}
				taxaList.add(new Taxon(values[0]));
				for (int i = 0; i < ntraits; i++) {
					traitValues.get(i)[dataCount] = values[i + 1];
				}
				dataCount++;
			}
			
			lineCount++;
		}
		
		phenotypeReader.close();
		
		attributes.add(new TaxaAttribute(taxaList));
		for (int i = 0; i < ntraits; i++) attributes.add(new CategoricalAttribute(traitnames[i], traitValues.get(i)));
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private void processHeader(String headerLine, Map<String, ArrayList<String>> headerMap) {
		String[] header = headerLine.split("[<>=\\s]+");
		if ( header[1].toLowerCase().equals("header") && header[2].toLowerCase().equals("name")) {
			String name = header[3];
			ArrayList<String> values = new ArrayList<>();
			for (int i = 4; i < header.length; i++) {
				values.add(header[i]);
			}
			headerMap.put(name, values);
		} else throw new IllegalArgumentException("Improperly formatted Header: " + headerLine);
	}
	
	private Phenotype processTraitsAndFactors(File phenotypeFile, String[] traitnames, int numberOfDataLines, boolean isCovariate, ArrayList<String> headerList)
	throws IOException {
		TreeSet<String> traitSet = new TreeSet<>();
		for (String trait : traitnames) traitSet.add(trait);
		HashMap<String, Integer> traitMap = new HashMap<>();
		int traitCount = 0;
		for (String trait : traitSet) traitMap.put(trait, traitCount++);
		
		int ntraitnames = traitnames.length;
		int ntraits = traitSet.size();
		int nfactors = headerList.size();
		String[] factorNames = new String[nfactors];
		
		//create set of composite headers and get factor names
		String[] splitHeader = headerList.get(0).split("[<>=\\s]+");
		factorNames[0] = splitHeader[3];
		String[] factorValues = Arrays.copyOfRange(splitHeader, 4, splitHeader.length);
		
		for (int i = 1; i < nfactors; i++) {
			splitHeader = headerList.get(i).split("[<>=\\s]+");
			factorNames[i] = splitHeader[3].replace("|", ":");
			for (int j = 0; j < factorValues.length; j++) {
				factorValues[j] += "|" + splitHeader[j + 4].replace("|", ":");
			}
		}
		
		TreeSet<String> factorSet = new TreeSet<>();
		for (String val:factorValues) factorSet.add(val);
		int nCompositeFactors = factorSet.size();
		HashMap<String,Integer> factorMap = new HashMap<>();
		int factorCount = 0;
		for (String factor : factorSet) factorMap.put(factor, factorCount++);
		
		//the length of each attribute array will be numberOfDataLines * number of composite headers
		int ndata = numberOfDataLines * nCompositeFactors;

		//create the factor arrays
		ArrayList<String[]> factorAttributeArrays = new ArrayList<>(nfactors);
		for (int i = 0; i < nfactors; i++) factorAttributeArrays.add(new String[ndata]);
		int fromIndex = 0;
		int toIndex = numberOfDataLines;
		for (String factor : factorSet) {
			String[] subFactor = factor.split("|");
			for (int i = 0; i < nfactors; i++) {
				Arrays.fill(factorAttributeArrays.get(i), fromIndex, toIndex, subFactor[i]);
			}
			fromIndex += numberOfDataLines;
			toIndex += numberOfDataLines;
		}
		
		//create the trait arrays and temporary taxa list
		ArrayList<float[]> traitAttributeArrays = new ArrayList<>();
		ArrayList<BitSet> missingList = new ArrayList<>();
		ArrayList<Taxon> tempTaxa = new ArrayList<>();

		for (int i = 0; i < ntraits; i++) {  //initialize the trait arrays to NaN
			float[] traitArray = new float[ndata];
			Arrays.fill(traitArray, Float.NaN);
			traitAttributeArrays.add(traitArray);
			missingList.add(new OpenBitSet(ndata));
		}
		
		BufferedReader br = new BufferedReader(new FileReader(phenotypeFile));
		int dataCount = 0;
		int lineCount = 1;
		String inputline;
		while ((inputline = br.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				String[] values = inputline.split("\\s+");
				if (values.length != ntraitnames + 1) {
					String msg = String.format("Incorrect number of values in line %d of %s", lineCount, phenotypeFile.getName());
					br.close();
					throw new IllegalArgumentException(msg);
				}
				tempTaxa.add(new Taxon(values[0]));
				for (int i = 0; i < ntraitnames; i++) {
					float val;
					int traitnum = traitMap.get(traitnames[i]);
					int factornum = factorMap.get(factorValues[i]);
					int dataIndex = factornum * numberOfDataLines + dataCount;
					try {
						val = Float.parseFloat(values[i + 1]);
					} catch (NumberFormatException e) {
						val = Float.NaN;
						missingList.get(traitnum).fastSet(dataIndex);
					}
					traitAttributeArrays.get(traitnum)[dataIndex] = val;
				}
				dataCount++;
			}
			
			lineCount++;
		}

		br.close();
		
		//create the taxa list
		ArrayList<Taxon> taxaList = new ArrayList<>(ndata);
		for (int i = 0; i < nCompositeFactors; i++) taxaList.addAll(tempTaxa);

		//create the taxa, factor, and trait attributes and types in that order
		ATTRIBUTE_TYPE myAttributeType;
		if (isCovariate) myAttributeType = ATTRIBUTE_TYPE.covariate;
		else myAttributeType = ATTRIBUTE_TYPE.data;

		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(nfactors + ntraits);
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(nfactors + ntraits);
		
		attributes.add(new TaxaAttribute(taxaList));
		types.add(ATTRIBUTE_TYPE.taxa);
		
		for (int i = 0; i < nfactors; i++) {
			attributes.add(new CategoricalAttribute(factorNames[i], factorAttributeArrays.get(i)));
			types.add(ATTRIBUTE_TYPE.factor);
		}
		
		traitCount = 0;
		for (String trait : traitSet) {
			attributes.add(new NumericAttribute(trait, traitAttributeArrays.get(traitCount), missingList.get(traitCount)));
			types.add(myAttributeType);
			traitCount++;
		}
		
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private Phenotype filterBasePhenotype() {
		if (attributeList != null) {
			if (attributeTypeList == null) {
				attributeTypeList = new ArrayList<ATTRIBUTE_TYPE>();
				for (PhenotypeAttribute attr:attributeList) {
					attributeTypeList.add(basePhenotype.attributeType(basePhenotype.indexOfAttribute(attr)));
				}
			}
			applyAttributeChangeMap();
			basePhenotype = new CorePhenotype(attributeList, attributeTypeList, phenotypeName);
			
		} else if (indexOfAttributesToKeep != null) {
			attributeList = new ArrayList<PhenotypeAttribute>();
			attributeTypeList = new ArrayList<ATTRIBUTE_TYPE>();
			for (int attrnum : indexOfAttributesToKeep) {
				attributeList.add(basePhenotype.attribute(attrnum));
			}
			if (attributeTypeList == null || attributeTypeList.size() != attributeList.size()) {
				for (int attrnum : indexOfAttributesToKeep) {
					attributeTypeList.add(basePhenotype.attributeType(attrnum));
				}
			}
			applyAttributeChangeMap();
			basePhenotype = new CorePhenotype(attributeList, attributeTypeList, phenotypeName);
			
		} else if (attributeChangeMap != null) {
			attributeList = basePhenotype.attributeListCopy();
			attributeTypeList = basePhenotype.typeListCopy();
			applyAttributeChangeMap();
			basePhenotype = new CorePhenotype(attributeList, attributeTypeList, phenotypeName);
		}

		if (taxaToKeep != null) {
			return FilterPhenotype.getInstance(basePhenotype, taxaToKeep, phenotypeName);
		}
		
		if (taxaToRemove != null) {
			TaxaList myTaxaList = basePhenotype.taxa();
			Iterator<Taxon> taxaIter = myTaxaList.iterator();
			while (taxaIter.hasNext()) {
				if (taxaToRemove.contains(taxaIter.next())) taxaIter.remove();
			}
			return FilterPhenotype.getInstance(basePhenotype, myTaxaList, phenotypeName);
		}

		return basePhenotype;
	}
	
	private void applyAttributeChangeMap() {
		if (attributeChangeMap != null && attributeTypeList != null && attributeList != null) {
			for (PhenotypeAttribute attr : attributeChangeMap.keySet()) {
				int ndx = attributeList.indexOf(attr);
				if (ndx > -1) attributeTypeList.set(ndx, attributeChangeMap.get(attr));
			}
		}
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
		if (phenotypesToJoin.size() < 2) throw new IllegalArgumentException("No join will be made because joining phenotypes requires at least two phenotypes.");
		if (isConcatenate) return concatenatePhenotypes();
		
		Iterator<Phenotype> phenoIter = phenotypesToJoin.iterator();
		Phenotype firstPhenotype = phenoIter.next();
		Phenotype secondPhenotype = phenoIter.next();
		Phenotype mergedPhenotype = mergeTwoPhenotypes(firstPhenotype, secondPhenotype);
		while (phenoIter.hasNext()) mergedPhenotype = mergeTwoPhenotypes(mergedPhenotype, phenoIter.next());
		return mergedPhenotype;
	}
	
	private Phenotype concatenatePhenotypes() {
		int nPheno = phenotypesToJoin.size();
		if (nPheno < 2) throw new IllegalArgumentException("No phenotypes to join.");
		attributeList = new ArrayList<PhenotypeAttribute>();
		attributeTypeList = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		
		//create a new TaxaAttribute
		List<Taxon> jointTaxaList = new ArrayList<Taxon>();
		String attrName = phenotypesToJoin.get(0).name();
		for (Phenotype pheno : phenotypesToJoin) {
			TaxaAttribute someTaxa = pheno.taxaAttribute();
			if (someTaxa == null) throw new IllegalArgumentException(String.format("Phenotypes cannot be concatenated because %s has no taxa.", pheno.name()));
			jointTaxaList.addAll(someTaxa.allTaxaAsList());
		}
		TaxaAttribute myTaxaAttribute = new TaxaAttribute(jointTaxaList, attrName);
		attributeList.add(myTaxaAttribute);
		attributeTypeList.add(ATTRIBUTE_TYPE.taxa);
		
		//make a list of attribute names (except for taxa) and find out how many total observations there will be
		HashSet<String> attributeNameSet = new HashSet<>();
		int totalNumberOfObs = 0;
		for (Phenotype pheno : phenotypesToJoin) {
			int nattr = pheno.numberOfAttributes();
			totalNumberOfObs += pheno.numberOfObservations();
			for (int i = 0; i < nattr; i++) {
				if (pheno.attributeType(i) != ATTRIBUTE_TYPE.taxa) attributeNameSet.add(pheno.attributeName(i));
			}
		}
		
		//create a new attribute for each name. Throw an error if attributes of the same name are of incompatible types.
		for (String attributeName : attributeNameSet) {
			ATTRIBUTE_TYPE myType = null;
			//take the type from the first phenotype with this name
			for (Phenotype pheno : phenotypesToJoin) {
				int attrIndex = pheno.attributeIndexForName(attributeName);
				if (attrIndex > -1) {
					myType = pheno.attributeType(attrIndex);
					break;
				}
			}
			
			//create this attribute
			if (myType == ATTRIBUTE_TYPE.factor) {
				String[] values = new String[totalNumberOfObs];
				int nPreviousObs = 0;
				for (Phenotype pheno : phenotypesToJoin) {
					int nObs = pheno.numberOfObservations();
					int attrIndex = pheno.attributeIndexForName(attributeName);
					if (attrIndex > -1) {
						System.arraycopy((String[]) pheno.attribute(attrIndex).allValues(), 0, values, nPreviousObs, nObs);
					} else {
						Arrays.fill(values, nPreviousObs, nPreviousObs + nObs, CategoricalAttribute.missingValue);
					}
					nPreviousObs += nObs;
				}
				attributeList.add(new CategoricalAttribute(attributeName, values));
				attributeTypeList.add(myType);
			} else {
				float[] values = new float[totalNumberOfObs];
				BitSet missing = new OpenBitSet(totalNumberOfObs);
				int nPreviousObs = 0;
				for (Phenotype pheno : phenotypesToJoin) {
					int nObs = pheno.numberOfObservations();
					int attrIndex = pheno.attributeIndexForName(attributeName);
					if (attrIndex > -1) {
						PhenotypeAttribute thisAttribute = pheno.attribute(attrIndex);
						System.arraycopy((float[]) thisAttribute.allValues(), 0, values, nPreviousObs, nObs);
						for (int i = 0; i < nObs; i++) if (thisAttribute.isMissing(i)) missing.fastSet(nPreviousObs + i);
					} else {
						Arrays.fill(values, nPreviousObs, nPreviousObs + nObs, Float.NaN);
						missing.set(nPreviousObs, nPreviousObs + nObs);
					}
					nPreviousObs += nObs;
				}
				attributeList.add(new NumericAttribute(attributeName, values, missing));
				attributeTypeList.add(myType);
			}
		}
		
		sortAttributes();
		
		return new CorePhenotype(attributeList, attributeTypeList, phenotypeName);
	}
	
	private Phenotype mergeTwoPhenotypes(Phenotype pheno1, Phenotype pheno2) {
		
		//build attribute name list for the new phenotype
		TreeSet<String> attributeNameSet = new TreeSet<>();
		
		int nAttributes1 = pheno1.numberOfAttributes();
		for (int a = 0; a < nAttributes1; a++) {
			if (pheno1.attributeType(a) != ATTRIBUTE_TYPE.taxa) attributeNameSet.add(pheno1.attributeName(a));
		}

		int nAttributes2 = pheno2.numberOfAttributes();
		for (int a = 0; a < nAttributes2; a++) {
			if (pheno2.attributeType(a) != ATTRIBUTE_TYPE.taxa) attributeNameSet.add(pheno2.attributeName(a));
		}
		
		ArrayList<String> attributeNameList = new ArrayList<>(attributeNameSet);
		ArrayList<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
		
		//create a type list for these attribute names
		for (String name : attributeNameList) {
			 int ndx = pheno1.attributeIndexForName(name);
			 if (ndx > -1) typeList.add(pheno1.attributeType(ndx));
			 else {
				 typeList.add(pheno2.attributeType(pheno2.attributeIndexForName(name)));
			 }
		}
		
		//create a list of taxa to be included in the merged phenotype
		TaxaList outTaxa;
		if (isUnionJoin) outTaxa = TaxaListUtils.getAllTaxa(pheno1.taxa(), pheno2.taxa());
		else outTaxa = TaxaListUtils.getCommonTaxa(pheno1.taxa(), pheno2.taxa());
		
		//for each phenotype create a Multimap with Taxon as key, obs number as value
		Multimap<Taxon, Integer> pheno1ObservationMap = HashMultimap.create();
		List<Taxon> pheno1Taxa = pheno1.taxaAttribute().allTaxaAsList();
		int taxonCount = 0;
		for (Taxon taxon:pheno1Taxa) pheno1ObservationMap.put(taxon, taxonCount++);
		
		Multimap<Taxon, Integer> pheno2ObservationMap = HashMultimap.create();
		List<Taxon> pheno2Taxa = pheno2.taxaAttribute().allTaxaAsList();
		taxonCount = 0;
		for (Taxon taxon:pheno2Taxa) pheno2ObservationMap.put(taxon, taxonCount++);
		
		//create a list of new observations
		//first find the factors (if any) in common
		//for each factor in common record index in pheno1 and pheno2
		ArrayList<String> mergeFactors = new ArrayList<String>();
		ArrayList<int[]> mergeFactorIndex = new ArrayList<>();
		int numberOfMergeFactors;
		if (pheno1.numberOfAttributesOfType(ATTRIBUTE_TYPE.factor) == 0 || pheno2.numberOfAttributesOfType(ATTRIBUTE_TYPE.factor) == 0) numberOfMergeFactors = 0;
		else {
			int[] factorIndex = pheno1.attributeIndicesOfType(ATTRIBUTE_TYPE.factor);
			for (int ndx1:factorIndex) {
				String factorName = pheno1.attributeName(ndx1);
				int ndx2 = pheno2.attributeIndexForName(factorName);
				if (ndx2 > -1 && pheno2.attributeType(ndx2) == ATTRIBUTE_TYPE.factor) {
					mergeFactors.add(factorName);
					mergeFactorIndex.add(new int[]{ndx1, ndx2});
				}
			}
		}
		numberOfMergeFactors = mergeFactors.size();
		
		if (numberOfMergeFactors > 0) {
			Collections.sort(mergeFactors);
			for (String factorName : mergeFactors) {
				mergeFactorIndex.add(new int[]{pheno1.attributeIndexForName(factorName), pheno2.attributeIndexForName(factorName)});
			}
		}
		
		//code for creating a list of new observations
		ArrayList<int[]> mergeObservation = new ArrayList<>();
		ArrayList<Taxon> listOfTaxaForObservations = new ArrayList<>();
		for (Taxon taxon:outTaxa) {
			Collection<Integer> pheno1obs = pheno1ObservationMap.get(taxon);
			Collection<Integer> pheno2obs = pheno2ObservationMap.get(taxon);
			boolean[] wasPheno2ObsUsed = new boolean[pheno2obs.size()];
			Arrays.fill(wasPheno2ObsUsed, false);
			for (Integer obs1 : pheno1obs) {
				boolean wasPheno1ObsUsed = false;
				int obs2Count = 0;
				for (Integer obs2 : pheno2obs) {
					boolean mergeTheseObs = true;
					for (int[] ndx : mergeFactorIndex) {
						if (pheno1.value(obs1, ndx[0]).equals(pheno2.value(obs1, ndx[1]))) mergeTheseObs = false;
						break;
					}
					if (mergeTheseObs) {
						wasPheno1ObsUsed = true;
						wasPheno2ObsUsed[obs2Count] = true;
						mergeObservation.add(new int[]{obs1, obs2});
						listOfTaxaForObservations.add(taxon);
					}
					obs2Count++;
				}
				if (!wasPheno1ObsUsed) {
					mergeObservation.add(new int[]{obs1, -1});
					listOfTaxaForObservations.add(taxon);
				}
			}
			int obs2Count = 0;
			for (Integer obs2:pheno2obs) {
				if (!wasPheno2ObsUsed[obs2Count++]) {
					mergeObservation.add(new int[]{-1, obs2});
					listOfTaxaForObservations.add(taxon);
				}
			}
		}
		
		//create the new attribute and type lists. Add the taxa attribute to them.
		TaxaAttribute newTaxaAttribute = new TaxaAttribute(listOfTaxaForObservations, pheno1.taxaAttribute().name());
		ArrayList<PhenotypeAttribute> newAttributes = new ArrayList<>();
		ArrayList<ATTRIBUTE_TYPE> newTypes = new ArrayList<>();
		newAttributes.add(newTaxaAttribute);
		newTypes.add(ATTRIBUTE_TYPE.taxa);
		
		//create and add the other attributes (in the name list) to the new attribute and type lists
		int nAdd = attributeNameList.size();
		int nObs = mergeObservation.size();
		for (int i = 0; i < nAdd; i++) {
			ATTRIBUTE_TYPE myType = typeList.get(i);
			String attrName = attributeNameList.get(i);
			
			//the attrnum array stores the index of this attribute in pheno1 and pheno2
			int[] attrnum = new int[]{pheno1.attributeIndexForName(attrName), pheno2.attributeIndexForName(attrName)};
			
			switch(myType) {
			case data:
			case covariate:
				float[] myFloatData = new float[nObs];
				BitSet myMissing = new OpenBitSet(nObs);
				int obsCount = 0;
				for (int[] ndx : mergeObservation) {
					boolean pheno1HasNonmissingValue = attrnum[0] > -1 && ndx[0] > -1 & !pheno1.isMissing(ndx[0], attrnum[0]);
					boolean pheno2HasNonmissingValue = attrnum[1] > -1 && ndx[1] > -1 & !pheno1.isMissing(ndx[1], attrnum[1]);
					if (pheno1HasNonmissingValue && pheno2HasNonmissingValue) {
						throw new IllegalArgumentException("Data sets will not be joined because both phenotypes have values for " + attrName);
					} else if (pheno1HasNonmissingValue) {
						myFloatData[obsCount] = (Float) pheno1.value(ndx[0], attrnum[0]);
					} else if (pheno2HasNonmissingValue) {
						myFloatData[obsCount] = (Float) pheno2.value(ndx[1], attrnum[1]);
					} else {
						myFloatData[obsCount] = Float.NaN;
						myMissing.fastSet(obsCount);
					}
					obsCount++;
				}
				newAttributes.add(new NumericAttribute(attrName, myFloatData, myMissing));
				newTypes.add(myType);
				break;
			case factor:
				boolean isMergeAttribute;
				if (mergeFactors.contains(attrName)) isMergeAttribute = true;
				else isMergeAttribute = false;

				String[] myStringData = new String[nObs];
				obsCount = 0;
				for (int[] ndx : mergeObservation) {
					boolean pheno1HasNonmissingValue = attrnum[0] > -1 && ndx[0] > -1 & !pheno1.isMissing(ndx[0], attrnum[0]);
					boolean pheno2HasNonmissingValue = attrnum[1] > -1 && ndx[1] > -1 & !pheno1.isMissing(ndx[1], attrnum[1]);
					if (pheno1HasNonmissingValue && pheno2HasNonmissingValue) {
						if (isMergeAttribute) myStringData[obsCount] = (String) pheno1.value(ndx[0], attrnum[0]);
						else throw new IllegalArgumentException("Data sets will not be joined because both phenotypes have values for " + attrName);
					} else if (pheno1HasNonmissingValue) {
						myStringData[obsCount] = (String) pheno1.value(ndx[0], attrnum[0]);
					} else if (pheno2HasNonmissingValue) {
						myStringData[obsCount] = (String) pheno2.value(ndx[1], attrnum[1]);
					} else {
						myStringData[obsCount] = CategoricalAttribute.missingValue;
					}
					obsCount++;
				}
				newAttributes.add(new CategoricalAttribute(attrName, myStringData));
				newTypes.add(myType);
				break;
			}
		}
		
		return new CorePhenotype(newAttributes, newTypes, phenotypeName);
	}
	
	private void sortAttributes() {
		//sort attributes, taxa first, then factors, then covariates, then data. Sorted by name within type
		int nAttr = attributeList.size();
		ArrayList<Integer> index = new ArrayList<>();
		for (int i = 0; i < nAttr; i++) index.add(i);
		
		Collections.sort(index, new Comparator<Integer>() {

			@Override
			public int compare(Integer o1, Integer o2) {
				ATTRIBUTE_TYPE type1 = attributeTypeList.get(o1);
				ATTRIBUTE_TYPE type2 = attributeTypeList.get(o2);
				if (type1 == type2) {
					return attributeList.get(o1).name().compareTo(attributeList.get(o2).name());
				} else {
					if (type1 == ATTRIBUTE_TYPE.taxa) return -1;
					else if (type2 == ATTRIBUTE_TYPE.taxa) return 1;
					else if (type1 == ATTRIBUTE_TYPE.factor) return -1;
					else if (type2 == ATTRIBUTE_TYPE.factor) return 1;
					else if (type1 == ATTRIBUTE_TYPE.covariate) return -1;
					else if (type2 == ATTRIBUTE_TYPE.covariate) return 1;
				}
				return 0;
			}
			
		});
		
		List<PhenotypeAttribute> unsortedAttributeList = attributeList;
		List<ATTRIBUTE_TYPE> unsortedAttributeTypeList = attributeTypeList;
		attributeList = new ArrayList<PhenotypeAttribute>();
		attributeTypeList = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		for (Integer ndx : index) {
			attributeList.add(unsortedAttributeList.get(ndx));
			attributeTypeList.add(unsortedAttributeTypeList.get(ndx));
		}
		
	}
	
}
