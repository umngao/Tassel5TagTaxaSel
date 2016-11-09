package net.maizegenetics.analysis.numericaltransform;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.taxa.Taxon;

public class AvgPhenotype {

    public Phenotype averagePheno(Phenotype pheno) throws Exception {
        //Get the variable names and put it into a list
        ArrayList<PhenotypeAttribute> attributes = (ArrayList<PhenotypeAttribute>)pheno.attributeListOfType(ATTRIBUTE_TYPE.data);
        ArrayList<String> phenoNames = (ArrayList<String>)attributes.stream().map((phenoAttr) -> phenoAttr.name()).collect(Collectors.toList());
        
        
        //Loop through the List<Taxa> get the unique names
        ArrayList<Taxon> taxonList = (ArrayList<Taxon>)pheno.taxaAttribute().allTaxaAsList();
        ArrayList<String> uniqueTaxaNames =(ArrayList<String>) pheno.taxaAttribute().allTaxaAsList().stream()
                                                                .map((taxon) -> taxon.getName()).distinct().collect(Collectors.toList());
        
        HashMap<String, ArrayList<ArrayList<Double>>> rowMapByTaxa= new HashMap<String,ArrayList<ArrayList<Double>>>();
        
        //Create an ArrayList object for each unique taxaName
        for(String taxaName : uniqueTaxaNames) {
            ArrayList<ArrayList<Double>> list = new ArrayList<ArrayList<Double>>();
            rowMapByTaxa.put(taxaName, list);
        }
        //Go through all taxon and use the name to associate duplicates together
        for(int i = 0; i < taxonList.size(); i++) {
            //Make a Container to hold the values
            ArrayList<Double> currentRow = new ArrayList<Double>();
            for(String varName : phenoNames) {
                //Add this to the duplicate list
                currentRow.add(((Float)pheno.value(i,pheno.attributeIndexForName(varName))).doubleValue());
            }
            //Add the duplicate row(by taxa name) to the mapping
            rowMapByTaxa.get(taxonList.get(i)).add(currentRow);
        }
        
        //Create 2-D ArrayList to store all of the averaged values
        ArrayList<ArrayList<Double>> avgValues = new ArrayList<ArrayList<Double>>();
        //for each taxa sum up all values and divide
        //Use Commons Math in case we want to do a geometric mean
        for(String taxaName : uniqueTaxaNames) {
            //Get the list of all the duplicate values for the taxaName
            ArrayList<ArrayList<Double>> currentValues = rowMapByTaxa.get(taxaName);
            
            //Make the container to hold the current average values
            ArrayList<Double> currentAverages = new ArrayList<Double>();
            //Loop through this 2-D array to and apply the average
            for(int col = 0; col < currentValues.get(0).size(); col++) {
                double[] values = new double[currentValues.size()];
                Mean meanObj = new Mean();
                for(int row = 0; row < currentValues.size(); row++) {
                   values[row] = currentValues.get(row).get(col);
                }
                //Add the averaged value to the finished current row
                currentAverages.add(meanObj.evaluate(values));
            }
            //Add the averaged current row to the avgValues 2-D array
            avgValues.add(currentAverages);
        }
        
        try {
            //Take the taxa List, the variable name list and the matrix and send to PhenoUtils package to get a phenotype back
            return PhenotypeUtils.createPhenotypeFromTransform(uniqueTaxaNames, phenoNames, avgValues);
        }
        catch(Exception e) {
            throw e;
        }
    }
}
