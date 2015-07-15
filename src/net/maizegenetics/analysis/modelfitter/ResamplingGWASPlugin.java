package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import javax.swing.ImageIcon;

import net.maizegenetics.analysis.modelfitter.ResidualForwardRegression;
import net.maizegenetics.analysis.modelfitter.AdditiveModelForwardRegression;
import net.maizegenetics.analysis.modelfitter.ForwardRegression;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.linearmodels.BasicShuffler;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public class ResamplingGWASPlugin extends AbstractPlugin {
    //input parameters :
    //enterlimit = the p-value cutoff for a term to enter the model
    //select next term using residuals
    //number of iterations
    //resample method (with/without replacement)
    //resample proportion
    //maximum number of terms in the model
    
    //restrictions:
    //there must be no missing data in the Phenotype
    //the input data set must contain one GenotypeTable and one Phenotype
    //at most one factor is allowed (expected to be family or pop)

    private Random randomGen = new Random();
    
    private PluginParameter<Double> enterLimit = new PluginParameter.Builder<>("enterLimit", 1e-8, Double.class)
            .description("A new term entering the model must have a p-value equal to or less than the enter limit. (Default = 1e-8)")
            .guiName("Enter Limit")
            .build();
    private PluginParameter<Integer> maxModelTerms = new PluginParameter.Builder<>("maxTerms", 100, Integer.class)
            .description("The maximum number of variants that will be fit. If the chromosome residuals are being fit, the maximum number of variants fit per chromosome. (Default = 100)")
            .guiName("Max terms")
            .build();
    
    private PluginParameter<Boolean> useResiduals = new PluginParameter.Builder<>("residuals", true, Boolean.class)
            .description("Should new terms be tested using residuals from the prior model? The analysis runs faster using this option. (Default = true)")
            .guiName("Use residuals")
            .build();
    
    private PluginParameter<Integer> numberOfIterations = new PluginParameter.Builder<>("nIterations", 100, Integer.class)
            .description("The number of times the data should be resampled. (Default = 100)")
            .guiName("Number of Iterations")
            .build();
    
    private PluginParameter<Double> resampleProportion = new PluginParameter.Builder<>("resampleProportion", 0.8, Double.class)
            .description("The size of the sample is resample proportion times the number of observations in the complete data. For bootstrap, set this value to 1 and with replacement to true. (Default = 0.8)")
            .guiName("")
            .build();
    
    private PluginParameter<Boolean> withReplacement = new PluginParameter.Builder<>("replacement", false, Boolean.class)
            .description("Should the sample be formed by sampling with replacement?  For bootstrap, set resample proportion to 1 and this value to true. (Default = false)")
            .guiName("With Replacement")
            .build();
    
    public ResamplingGWASPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    
    @Override
    public DataSet processData(DataSet input) {
        //check that there is only one GenotypePhenotype
        List<Datum> datumList = input.getDataOfType(GenotypePhenotype.class);
        if (datumList.size() != 1) {
            throw new IllegalArgumentException("Exactly one joined Genotype-Phenotype dataset must be supplied as input to Resample GWAS.");
        }
        
        GenotypePhenotype myGenoPheno = (GenotypePhenotype) datumList.get(0).getData();
        
        //if the Phenotype has more than one factor level, throw an error
        int numberOfFactors = myGenoPheno.phenotype().numberOfAttributesOfType(ATTRIBUTE_TYPE.factor);
        if (numberOfFactors > 1) 
            throw new IllegalArgumentException("Phenotype supplied to Resample GWAS can have at most one factor");
        
        //if there is any missing data in the phenotype throw an error
        boolean anyMissing = myGenoPheno.phenotype().attributeStream().anyMatch(pa -> pa.missing().cardinality() > 0 );
        if (anyMissing) 
            throw new IllegalArgumentException("No missing phenotype data allow as input to Resample GWAS");
        
        int numberOfTraits = myGenoPheno.phenotype().numberOfAttributesOfType(ATTRIBUTE_TYPE.data);
        List<int[]> factorLevelList = new ArrayList<>();
        
        //create a factor level list for creating subsamples
        if (numberOfFactors == 1) {
            CategoricalAttribute myFactor = (CategoricalAttribute) myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor).get(0);
            int[] levels = myFactor.allIntValues();
            factorLevelList = new ArrayList<>();
            Map<Integer, List<Integer>> factorGroups =  IntStream.range(0, levels.length).boxed().collect(Collectors.groupingBy(i -> new Integer(levels[i])));
            for (int i = 0; i < factorGroups.size(); i++) {
                List<Integer> factorMembers = factorGroups.get(i);
                factorLevelList.add( factorMembers.stream().mapToInt(I -> I.intValue()).toArray() );
            }
        } else {
            factorLevelList.add(IntStream.range(0,  myGenoPheno.phenotype().numberOfObservations()).toArray());
        }
        
        String dataname = datumList.get(0).getName();
        TableReportBuilder reportBuilder = TableReportBuilder.getInstance("Resample terms_" + dataname, ForwardRegression.columnLabels());
        
        //for each chromosome
        for (int ph = 0; ph < numberOfTraits; ph++) {
            ForwardRegression modelfitter;
            if (useResiduals.value()) {
                modelfitter = new ResidualForwardRegression(myGenoPheno, ph, enterLimit.value(), maxModelTerms.value());
            } else {
                modelfitter = new AdditiveModelForwardRegression(myGenoPheno, ph, enterLimit.value(), maxModelTerms.value());
            }

            //for each iteration
            for (int iter = 0; iter < numberOfIterations.value(); iter++) {
                //create a random subsample
                int[] subsample = randomSample(factorLevelList);
                
                //generate a model
                modelfitter.fitModelForSubsample(subsample);
                
                //add the model terms and p-values to the result TableReportBuilder
                for (Object[] row : modelfitter.fittedModel()) reportBuilder.add(row);
            }
            
        }
        
        //return a DataSet that contains the TableReportBuilder
        Datum theResult = new Datum("name", reportBuilder.build(), "comment");
        return new DataSet(theResult, this);
    }

    private int[] randomSample(List<int[]> factorLevelList) {
        //if there is a single factor, will draw samples within factor levels
        double rp = resampleProportion.value().doubleValue();
        if (withReplacement.value()) {
            return factorLevelList.stream()
                    .flatMapToInt(iarray -> randomGen.ints((int) Math.round(iarray.length * rp), 0, iarray.length).map(i -> iarray[i]))
                    .toArray();
        } else {
            return factorLevelList.stream()
                    .flatMapToInt(iarray -> {
                        int[] copy = Arrays.copyOf(iarray, iarray.length);
                        BasicShuffler.shuffle(copy);
                        return Arrays.stream(copy).limit((int) Math.round(iarray.length * rp));
                    })
                    .toArray();
        }
    }

    @Override
    public String pluginDescription() {
        return "ResamplingGWASPlugin uses forward regression to fit a multiple SNP model to each of a number of samples drawn from the phenotype data. ";
    }


    @Override
    public ImageIcon getIcon() {
        URL imageURL = FixedEffectLMPlugin.class.getResource("/net/maizegenetics/analysis/images/resample.png");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Resample GWAS";
    }

    @Override
    public String getToolTipText() {
        return "GWAS with resampling";
    }
    
    /**
     * Initializes the random number generator with a seed. Primarily used for testing to generate reproducible results.
     * @param seed      the seed for the Random object
     */
    public void setRandomSeed(int seed) {
        randomGen = new Random(seed);
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ResamplingGWASPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

    /**
     * A new term entering the model must have a p-value equal
     * to or less than the enter limit. (Default = 1e-8)
     *
     * @return Enter Limit
     */
    public Double enterLimit() {
        return enterLimit.value();
    }

    /**
     * Set Enter Limit. A new term entering the model must
     * have a p-value equal to or less than the enter limit.
     * (Default = 1e-8)
     *
     * @param value Enter Limit
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin enterLimit(Double value) {
        enterLimit = new PluginParameter<>(enterLimit, value);
        return this;
    }

    /**
     * Should new terms be tested using residuals from the
     * prior model? The analysis runs faster using this option.
     * (Default = true)
     *
     * @return Use residuals
     */
    public Boolean useResiduals() {
        return useResiduals.value();
    }

    /**
     * Set Use residuals. Should new terms be tested using
     * residuals from the prior model? The analysis runs faster
     * using this option. (Default = true)
     *
     * @param value Use residuals
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin useResiduals(Boolean value) {
        useResiduals = new PluginParameter<>(useResiduals, value);
        return this;
    }

    /**
     * The number of times the data should be resampled. (Default
     * = 100)
     *
     * @return Number of Iterations
     */
    public Integer numberOfIterations() {
        return numberOfIterations.value();
    }

    /**
     * Set Number of Iterations. The number of times the data
     * should be resampled. (Default = 100)
     *
     * @param value Number of Iterations
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin numberOfIterations(Integer value) {
        numberOfIterations = new PluginParameter<>(numberOfIterations, value);
        return this;
    }

    /**
     * The size of the sample is resample proportion times
     * the number of observations in the complete data. For
     * bootstrap, set this value to 1 and with replacement
     * to true. (Default = 0.8)
     *
     * @return Resample Proportion
     */
    public Double resampleProportion() {
        return resampleProportion.value();
    }

    /**
     * Set Resample Proportion. The size of the sample is
     * resample proportion times the number of observations
     * in the complete data. For bootstrap, set this value
     * to 1 and with replacement to true. (Default = 0.8)
     *
     * @param value Resample Proportion
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin resampleProportion(Double value) {
        resampleProportion = new PluginParameter<>(resampleProportion, value);
        return this;
    }

    /**
     * Should the sample be formed by sampling with replacement?
     * (Default = false)
     *
     * @return With Replacement
     */
    public Boolean withReplacement() {
        return withReplacement.value();
    }

    /**
     * Set With Replacement. Should the sample be formed by
     * sampling with replacement? (Default = false)
     *
     * @param value With Replacement
     *
     * @return this plugin
     */
    public ResamplingGWASPlugin withReplacement(Boolean value) {
        withReplacement = new PluginParameter<>(withReplacement, value);
        return this;
    }

}
