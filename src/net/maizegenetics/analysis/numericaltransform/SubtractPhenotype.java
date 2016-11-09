package net.maizegenetics.analysis.numericaltransform;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Collectors;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;

public class SubtractPhenotype {
    //Method to subtract phenotypes where the taxa Names correspond
    public Phenotype subtractPhenotype(Phenotype pheno1, Phenotype pheno2) throws Exception {
        //TODO verify no duplicate taxa otherwise it will take the last taxa in the list for subtraction
        
        //Get the intersection of taxa names between pheno1 and pheno2
        ArrayList<String> uniqueTaxaPheno1 = (ArrayList<String>)pheno1.taxaAttribute().allTaxaAsList().stream().map((taxon) -> taxon.getName())
                                                                        .distinct().collect(Collectors.toList());
        ArrayList<String> uniqueTaxaPheno2 = (ArrayList<String>)pheno2.taxaAttribute().allTaxaAsList().stream().map((taxon) -> taxon.getName())
                                                                        .distinct().collect(Collectors.toList());
        ArrayList<String> intersectionTaxaList = (ArrayList<String>) uniqueTaxaPheno1.clone();
        intersectionTaxaList.retainAll(uniqueTaxaPheno2);
        //Create mapping so we can get the row index for a given taxa
        HashMap<String,Integer> pheno1NameToIndexMap = new HashMap<String,Integer>();
        HashMap<String,Integer> pheno2NameToIndexMap = new HashMap<String,Integer>();

        //Loop through the taxa lists to get the index for each taxa
        for(int i = 0; i < uniqueTaxaPheno1.size(); i++) {
            pheno1NameToIndexMap.put(uniqueTaxaPheno1.get(i),i);
        }
        
        for(int i = 0; i < uniqueTaxaPheno2.size(); i++) {
            pheno2NameToIndexMap.put(uniqueTaxaPheno2.get(i),i);
        }
        
        //Get the intersection of Phenotype Variable names
        ArrayList<PhenotypeAttribute> attributes1 = (ArrayList<PhenotypeAttribute>)pheno1.attributeListOfType(ATTRIBUTE_TYPE.data);
        ArrayList<String> phenoNames1 = (ArrayList<String>)attributes1.stream().map((phenoAttr) -> phenoAttr.name()).collect(Collectors.toList());
        
        ArrayList<PhenotypeAttribute> attributes2 = (ArrayList<PhenotypeAttribute>)pheno2.attributeListOfType(ATTRIBUTE_TYPE.data);
        ArrayList<String> phenoNames2 = (ArrayList<String>)attributes2.stream().map((phenoAttr) -> phenoAttr.name()).collect(Collectors.toList());
        
        ArrayList<String> intersectionPhenoList = (ArrayList<String>) phenoNames1.clone();
        intersectionPhenoList.retainAll(phenoNames2);
        
        //Transform each phenotype object to be Map<taxaName, data>
        ArrayList<ArrayList<Double>> data = new ArrayList<ArrayList<Double>>();
        for(int row = 0; row < intersectionTaxaList.size(); row++) {
            ArrayList<Double> currentRow = new ArrayList<Double>();
            for(int col = 0; col < intersectionPhenoList.size(); col++) {
                //Get the value from the first phenotype object
                double firstVal = ((Float)pheno1.value(pheno1NameToIndexMap.get(intersectionTaxaList.get(row)), 
                            pheno1.attributeIndexForName(intersectionPhenoList.get(col)))).doubleValue();
                //Get the value from the second phenotype object
                double secondVal = ((Float)pheno2.value(pheno2NameToIndexMap.get(intersectionTaxaList.get(row)), 
                        pheno2.attributeIndexForName(intersectionPhenoList.get(col)))).doubleValue();
                currentRow.add(firstVal - secondVal);
            }
            data.add(currentRow);
        }
        try {
        //Take the taxa List, the variable name list and the matrix and send to PhenoUtils package to get a phenotype back
        return PhenotypeUtils.createPhenotypeFromTransform(intersectionTaxaList, intersectionPhenoList, data);
        }
        catch(Exception e) {
            throw e;
        }
    }
}
