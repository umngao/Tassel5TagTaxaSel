package net.maizegenetics.analysis.association;

import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import org.apache.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.EMMA.EMMAforDoubleMatrix;
import net.maizegenetics.stats.linearmodels.BasicShuffler;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class GenomicSelectionPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(RidgeRegressionEmmaPlugin.class);

    private PluginParameter<Boolean> performCrossValidation = new
        PluginParameter.Builder<>("doCV", true, Boolean.class)
            .description("Perform cross-validation: True or False")
            .guiName("Perform cross-validation")
            .build();
    
    private PluginParameter<Integer> kFolds = new
        PluginParameter.Builder<>("kFolds", 5, Integer.class)
            .description("Number of folds to use for k-fold cross-validation (default = 5)")
            .guiName("Number of folds")
            .build();
    
    private PluginParameter<Integer> nIterations = new
        PluginParameter.Builder<>("nIter", 20, Integer.class)
            .description("Number of iterations when running k-fold cross-validation (default = 20)")
            .guiName("Number of iterations")
            .build();
    
    public GenomicSelectionPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

//    @Override
//    public DataSet performFunction(DataSet input) {
//        try {
//            List<Datum> datasets = input.getDataOfType(GenotypePhenotype.class);
//            if (datasets.size() < 1) {
//                String msg = "No datasets of an appropriate type were selected for the GS analysis.";
//                myLogger.error(msg);
//                if (isInteractive()) {
//                    JOptionPane.showMessageDialog(getParentFrame(), msg, "GS Error", JOptionPane.ERROR_MESSAGE);
//                }
//                return null;
//            }
//
//            LinkedList<Datum> results = new LinkedList<Datum>();
//            for (Datum dataset : datasets) {
//                try {
//                    LinkedList<Datum> aResult = null;
//                    aResult = processData(dataset);
//                    if (aResult != null) {
//                        results.addAll(aResult);
//                        fireDataSetReturned(new DataSet(aResult, this));
//                    }
//                } catch (Exception e) {
//                    StringBuilder msg = new StringBuilder("Error in GS processing " + dataset.getName());
//                    msg.append(". ").append(e.getMessage());
//                    myLogger.error(msg.toString());
//                    e.printStackTrace();
//                    if (isInteractive()) {
//                        JOptionPane.showMessageDialog(getParentFrame(), msg.toString(), "GS Error", JOptionPane.ERROR_MESSAGE);
//                    }
//                }
//            }
//
//            return new DataSet(results, this);
//        } finally {
//            fireProgress(100);
//        }
//    }

    public DataSet processData(DataSet input) {
        DoubleMatrix phenotype;
        DistanceMatrix kinshipOriginal;
        DoubleMatrix fixedEffects;
        LinkedList<Datum> theResults = new LinkedList<Datum>();
        List<Datum> myDataList = input.getDataOfType(Phenotype.class);
        if (myDataList.size() == 0)
            throw new IllegalArgumentException("No phenotype selected.");
        if (myDataList.size() > 1)
            throw new IllegalArgumentException("Too many phenotypes selected.");
        Phenotype myPhenotype = (Phenotype) myDataList.get(0).getData();
        
        myDataList = input.getDataOfType(DistanceMatrix.class);
        if (myDataList.size() == 0)
            throw new IllegalArgumentException("No kinship matrix selected.");
        if (myDataList.size() > 1)
            throw new IllegalArgumentException("Too many kinship matrices selected.");
        kinshipOriginal = (DistanceMatrix) myDataList.get(0).getData();
        
          //Remove phenos for which no kinship is present
          TaxaList phenoTaxa = myPhenotype.taxa();
          TaxaList kinTaxa = kinshipOriginal.getTaxaList();
          TaxaList jointTaxa = TaxaListUtils.getCommonTaxa(phenoTaxa,kinTaxa);
          Phenotype reducedPheno = new PhenotypeBuilder()
                  .fromPhenotype(myPhenotype)
                  .keepTaxa(jointTaxa)
                  .build().get(0);
          
          //Remove kinship for which no pheno is present
          DistanceMatrix myKinship = DistanceMatrixUtils.keepTaxa(kinshipOriginal,jointTaxa);
          
        //numbers of different things
//        List<PhenotypeAttribute> dataAttributeList = reducedPheno.attributeListOfType(ATTRIBUTE_TYPE.data);
//        List<PhenotypeAttribute> factorAttributeList = reducedPheno.attributeListOfType(ATTRIBUTE_TYPE.factor);
//        List<PhenotypeAttribute> covariateAttributeList = reducedPheno.attributeListOfType(ATTRIBUTE_TYPE.covariate);
//        TaxaAttribute myTaxaAttribute = reducedPheno.taxaAttribute();

        //create an index that contains all attributes except the data attributes
        //the last index value will be unassigned and will take the value of each data attribute one at a time
        int[] taxaAttrIndex = reducedPheno.attributeIndicesOfType(ATTRIBUTE_TYPE.taxa);
        int[] dataAttrIndex = reducedPheno.attributeIndicesOfType(ATTRIBUTE_TYPE.data);
        int[] factorAttrIndex = reducedPheno.attributeIndicesOfType(ATTRIBUTE_TYPE.factor);
        int[] covariateAttrIndex = reducedPheno.attributeIndicesOfType(ATTRIBUTE_TYPE.covariate);
        int nAttributes = taxaAttrIndex.length + factorAttrIndex.length + covariateAttrIndex.length + 1;
        int[] singlePhenotypeIndex = new int[nAttributes];
        int counter = 0;
        for (int addIndex : taxaAttrIndex) singlePhenotypeIndex[counter++] = addIndex;
        for (int addIndex : factorAttrIndex) singlePhenotypeIndex[counter++] = addIndex;
        for (int addIndex : covariateAttrIndex) singlePhenotypeIndex[counter++] = addIndex;
        
        for (int singleDataIndex : dataAttrIndex) {
        	singlePhenotypeIndex[nAttributes - 1] = singleDataIndex;
        	Phenotype singlePhenotype = new PhenotypeBuilder().fromPhenotype(reducedPheno)
        			.keepAttributes(singlePhenotypeIndex).removeMissingObservations().build().get(0);
            List<PhenotypeAttribute> factorAttributeList = singlePhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);
            List<PhenotypeAttribute> covariateAttributeList = singlePhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);
//            TaxaAttribute myTaxaAttribute = singlePhenotype.taxaAttribute();

            
        	NumericAttribute dataAttribute = (NumericAttribute) singlePhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).get(0);
        	
            //get phenotype data
            double[] phenotypeData = AssociationUtils.convertFloatArrayToDouble(dataAttribute.floatValues());
            int nObs = phenotypeData.length;
            phenotype = DoubleMatrixFactory.DEFAULT.make(nObs, 1, phenotypeData);

            //make the fixed effect matrix
            int numberOfFactors = factorAttributeList.size();
            int numberOfCovariates = covariateAttributeList.size();
            int numberOfEffects = numberOfFactors + numberOfCovariates + 1;
            
            if (numberOfEffects > 1) {
                DoubleMatrix[][] effects = new DoubleMatrix[1][numberOfEffects];
                effects[0][0] = DoubleMatrixFactory.DEFAULT.make(nObs, 1, 1);
                for (int i = 0; i < numberOfFactors; i++) {
                	CategoricalAttribute fa = (CategoricalAttribute) factorAttributeList.get(i);
                    FactorModelEffect fme = new FactorModelEffect(fa.allIntValues(), true);
                    effects[0][i + 1] = fme.getX();
                }
                for (int i = 0; i < numberOfCovariates; i++) {
                	NumericAttribute na = (NumericAttribute) covariateAttributeList.get(i);
                	double[] values = AssociationUtils.convertFloatArrayToDouble(na.floatValues());
                    effects[0][i + numberOfFactors + 1] = DoubleMatrixFactory.DEFAULT.make(nObs, 1, values);
                }
                fixedEffects = DoubleMatrixFactory.DEFAULT.compose(effects);
            } else {
                fixedEffects = DoubleMatrixFactory.DEFAULT.make(nObs, 1, 1);
            }
            
            DoubleMatrix kinship = DoubleMatrixFactory.DEFAULT.make(myKinship.getClonedDistances());
            
            //Create folds for cross-validation
            //this division automatically rounds down. feature of java, tested to confirm.
            int foldSize = nObs/kFolds.value();
            int[] seq = IntStream.range(0,nObs).toArray();
            BasicShuffler.reset();//seed is 111 by default within BasicShuffler. This line returns to beginning of that sequence.
            int numberOfIterations = nIterations.value();
            int numberOfFolds = kFolds.value();
            double[] rValues = new double[numberOfIterations*numberOfFolds];
            int rValueIndex = 0;
            //Iterate through number of iterations specified.
            for (int iter = 0 ; iter < numberOfIterations ; iter++){
                BasicShuffler.shuffle(seq);
                //Iterate through folds for cross-validation
                int startFold = 0;
                for (int fold = 0 ; fold < numberOfFolds ; fold++){
                    //Set prediction fold to missing
                    DoubleMatrix phenoTraining = phenotype.copy();
                    int endFold = startFold + foldSize;
                    if (fold == numberOfFolds - 1) endFold = nObs;
                    for (int ndx = startFold ; ndx < endFold ; ndx++){
                        phenoTraining.set(seq[ndx],0,Double.NaN);
                    }
                    
                    //Run EMMA using G-BLUP constructor (data, fixed, kinship)
                    EMMAforDoubleMatrix runEMMA = new EMMAforDoubleMatrix(phenoTraining,fixedEffects,kinship);
                    runEMMA.solve();
                    runEMMA.calculateBlupsPredicted();
                    double[] predictions = runEMMA.getPred().to1DArray();
                    int testSize = endFold - startFold;
                    double[] testPredictions = new double[testSize];
                    double[] testObserved = new double[testSize];
                    for (int ndx = 0 ; ndx < testSize ; ndx++){
                        int seqIndex = seq[ndx + startFold];
                        testPredictions[ndx] = predictions[seqIndex];
                        testObserved[ndx] = phenotype.get(seqIndex,0);
                    }
                    PearsonsCorrelation Pearsons = new PearsonsCorrelation();
                    rValues[rValueIndex++] = Pearsons.correlation(testPredictions,testObserved);
                    startFold = endFold;
                }
            }
            
            DescriptiveStatistics stats = new DescriptiveStatistics(rValues);
            double meanR = stats.getMean();
            double varR = stats.getVariance();
            double sdR = Math.sqrt(varR/rValues.length);
            
            System.out.printf("For phenotype %s\n",dataAttribute.name());
            System.out.printf("Mean from genomic prediction = %1.4f\n",meanR);
            System.out.printf("Standard deviation of mean from genomic prediction = %1.8f\n",sdR);
            Arrays.stream(rValues).forEach(v -> System.out.printf("%1.4f\n",v));
            
            //StringBuilder comment = new StringBuilder("Genomic prediction with Ridge Regression");
            //comment.append("based on k-fold cross-validation with specified number of iterations.");
            //theResults.add(new Datum("Predictive Ability", predictiveAbility, comment.toString()));
            
//            float[] gebv;
//            gebv = AssociationUtils.convertDoubleArrayToFloat(predictions);
//            List<Integer> gebvIndex = IntStream.range(0,gebv.length).boxed().collect(Collectors.toList());
//            Collections.sort(gebvIndex, new Comparator<Integer>(){
//
//                    @Override
//                    public int compare(Integer t, Integer t1) {
//                        if (gebv[t] > gebv[t1]) return -1;
//                        if (gebv[t] < gebv[t1]) return 1;
//                        return 0;
//                    }      
//            });
//            
//            int[] sortedIndex = gebvIndex.stream().mapToInt(val -> val.intValue()).toArray();
//            
//            String phenoName = attr.name();
//            ArrayList<PhenotypeAttribute> attributeList = new ArrayList<>();
//            ArrayList<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
//            
//            PhenotypeAttribute myNewTaxa = myTaxaAttribute.subset(sortedIndex,myTaxaAttribute.name());
//            attributeList.add(myNewTaxa);
//            typeList.add(ATTRIBUTE_TYPE.taxa);
//            
//            PhenotypeAttribute myPheno = new NumericAttribute(phenoName + "_GEBV", gebv, new OpenBitSet(nObs));
//            myPheno = myPheno.subset(sortedIndex,myPheno.name());
//            attributeList.add(myPheno);
//            typeList.add(ATTRIBUTE_TYPE.data);
//            
//            Phenotype gebvPheno = new PhenotypeBuilder().assignName("GEBV_" + phenoName)
//            		.fromAttributeList(attributeList, typeList)
//            		.build().get(0);
//            
//            String datumName = "GEBVs_" + phenoName + "_"; //+ dataset.getName()
//            StringBuilder comment2 = new StringBuilder("Ridge Regression from ");
//            //comment2.append(dataset.getName()).append(":\n");
//            comment2.append("Genomic Estimated Breeding Values (GEBVs)\n");
//            comment2.append("trait = ").append(phenoName).append("\n");
//            comment2.append(nObs).append(" lines");
//            theResults.add(new Datum(datumName, gebvPheno, comment2.toString()));
        }
//
//        fireDataSetReturned(new DataSet(theResults, this));
//            	return new DataSet(theResults, this);
        return null;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = GenomicSelectionPlugin.class.getResource("/net/maizegenetics/analysis/images/LinearAssociation.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Genomic Selection";
    }

    @Override
    public String getToolTipText() {
        return "Predict Phenotypes using G-BLUP for Genomic Selection";
    }

    @Override
    public String pluginDescription() {
        return "Predicts phenotypes using G-BLUP for genomic selection using a user-inputted kinship matrix and phenotype(s).";
    }  
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GenomicSelectionPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     * @param input
     */
    // TODO: Replace <Type> with specific type.
    public DataSet runPlugin(DataSet input) {
        return (DataSet) performFunction(input).getData(0).getData();
    }

    /**
     * Perform cross-validation: True or False
     *
     * @return Perform cross-validation
     */
    public Boolean performCrossValidation() {
        return performCrossValidation.value();
    }

    /**
     * Set Perform cross-validation. Perform cross-validation:
     * True or False
     *
     * @param value Perform cross-validation
     *
     * @return this plugin
     */
    public GenomicSelectionPlugin performCrossValidation(Boolean value) {
        performCrossValidation = new PluginParameter<>(performCrossValidation, value);
        return this;
    }

    /**
     * Number of folds to use for k-fold cross-validation
     * (default = 5)
     *
     * @return Number of folds
     */
    public Integer kFolds() {
        return kFolds.value();
    }

    /**
     * Set Number of folds. Number of folds to use for k-fold
     * cross-validation (default = 5)
     *
     * @param value Number of folds
     *
     * @return this plugin
     */
    public GenomicSelectionPlugin kFolds(Integer value) {
        kFolds = new PluginParameter<>(kFolds, value);
        return this;
    }

    /**
     * Number of iterations when running k-fold cross-validation
     * (default = 20)
     *
     * @return Number of iterations
     */
    public Integer nIterations() {
        return nIterations.value();
    }

    /**
     * Set Number of iterations. Number of iterations when
     * running k-fold cross-validation (default = 20)
     *
     * @param value Number of iterations
     *
     * @return this plugin
     */
    public GenomicSelectionPlugin nIterations(Integer value) {
        nIterations = new PluginParameter<>(nIterations, value);
        return this;
    }
}
