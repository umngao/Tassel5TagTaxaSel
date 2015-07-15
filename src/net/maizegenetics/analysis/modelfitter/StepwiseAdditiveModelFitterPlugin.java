package net.maizegenetics.analysis.modelfitter;

import java.awt.Frame;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;

import com.google.common.collect.Range;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;

public class StepwiseAdditiveModelFitterPlugin extends AbstractPlugin {
    //needs to set the following parameters in StepwiseAdditiveModelFitter:
    //    private int numberOfPermutations = 1000;
    //    private double permutationAlpha = 0.05;
    //    private double enterLimit = 1e-5;
    //    private double exitLimit = 2e-5;
    //    private boolean useReferenceProbability = true;
    //    private boolean isNested = true;
    //    private String nestingEffectName = "family";
    //    private AdditiveSite.CRITERION modelSelectionCriterion = AdditiveSite.CRITERION.pval;
    //    private int maxSitesInModel = 1000;
    private List<String> myFactorNameList;
    private GenotypePhenotype myGenoPheno;
    private String datasetName;

    private PluginParameter<AdditiveSite.CRITERION> modelCriterion =
            new PluginParameter.Builder<>("criterion", AdditiveSite.CRITERION.pval, AdditiveSite.CRITERION.class)
                    .range(AdditiveSite.CRITERION.values())
                    .guiName("Model selection criterion")
                    .description("The model selection criterion used to determine which terms enter the model and how many. Value must be one of pval, bic, mbic, or aic")
                    .build();

    private PluginParameter<Boolean> useResiduals =
            new PluginParameter.Builder<>("useResidual", false, Boolean.class)
                    .description("Should the each new term be tested using residuals from the previous model instead of fitting a complete model each time?")
                    .guiName("Fit using residuals")
                    .build();

    private PluginParameter<Boolean> usePermutations =
            new PluginParameter.Builder<>("usePerm", true, Boolean.class)
                    .description("Should permutations be used to set the enter and exit limits for stepwise regression? A permutation test will be used to determine the enter limit. The exit limit will be set to 2 times the enter limit.")
                    .guiName("Use permutations")
                    .build();

    private PluginParameter<Integer> numberOfPermutations =
            new PluginParameter.Builder<>("nPerm", 1000, Integer.class)
                    .description("The number of permutations used to determine the enter limit.")
                    .guiName("Number of permutations")
                    .dependentOnParameter(usePermutations)
                    .build();

    private PluginParameter<Double> permutationAlpha =
            new PluginParameter.Builder<>("permAlpha", 0.05, Double.class)
                    .description("Type I errors will be controlled at this level.")
                    .guiName("Alpha for permutations")
                    .dependentOnParameter(usePermutations)
                    .build();

    private PluginParameter<Double> enterLimit =
            new PluginParameter.Builder<>("enterLimit", 1e-5, Double.class)
                    .description("When p-value is the model selection criteria, model fitting will stop when the next term chosen has a p-value greater than the enterLimit. This value will be over-ridden the permutation test, if used.")
                    .guiName("enterLimit")
                    .dependentOnParameter(usePermutations, false)
                    .build();

    private PluginParameter<Double> exitLimit =
            new PluginParameter.Builder<>("exitLimit", 2e-5, Double.class)
                    .description("During the backward step of model fitting if p-value has been chosen as the selection criterion, if the term in model with the highest p-value has a p-value > exitLimit, it will be removed from the model.")
                    .guiName("exitLimit")
                    .dependentOnParameter(usePermutations, false)
                    .build();

    private PluginParameter<Integer> maxTerms =
            new PluginParameter.Builder<>("maxSNPs", 1000, Integer.class)
                    .description("The maximum number of SNPs/markers that will be fit. If the model selection criterion is met first, then the fitting process will stop at that point.")
                    .guiName("Max SNPs/markers")
                    .build();

    private PluginParameter<Boolean> isNested =
            new PluginParameter.Builder<>("isNested", true, Boolean.class)
                    .description("Should SNPs/markers be nested within a factor, such as family?")
                    .guiName("")
                    .build();

    private PluginParameter<List> nestingFactor =
            new PluginParameter.Builder<>("nestFactor", null, List.class)
                    .guiName("Nesting factor")
                    .description("Nest markers within this factor.")
                    .dependentOnParameter(isNested)
                    .objectListSingleSelect()
                    .build();

    private PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myGenotypeTable =
            new PluginParameter.Builder<>("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.class)
                    .genotypeTable()
                    .range(GenotypeTable.GENOTYPE_TABLE_COMPONENT.values())
                    .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
                    .build();

    private PluginParameter<Boolean> createAnova =
            new PluginParameter.Builder<>("anova", true, Boolean.class)
                    .description("Create pre- and post-scan anova reports.")
                    .guiName("Create anova reports")
                    .build();

    private PluginParameter<Boolean> createEffects =
            new PluginParameter.Builder<>("effects", true, Boolean.class)
                    .description("Create a report of marker effects based on the scan results.")
                    .guiName("Create effects report")
                    .build();

    private PluginParameter<Boolean> createEffectsPrescan =
            new PluginParameter.Builder<>("effectsPrescan", false, Boolean.class)
                    .description("Create a report of marker effects based on the results pre-scan.")
                    .guiName("Create prescan effects")
                    .build();

    private PluginParameter<Boolean> createStep =
            new PluginParameter.Builder<>("step", true, Boolean.class)
                    .description("Create a report of the which markers enter and leave the model as it is being fit.")
                    .guiName("Create step report")
                    .build();

    private PluginParameter<Boolean> createResiduals =
            new PluginParameter.Builder<>("anova", false, Boolean.class)
                    .description("Create a phenotype dataset of model residuals for each chromosome. For each chromosome, the residuals will be calculated from a model with all terms EXCEPT the markers on that chromosome.")
                    .guiName("Create residuals")
                    .build();

    private PluginParameter<Boolean> writeFiles =
            new PluginParameter.Builder<>("saveToFile", false, Boolean.class)
                    .description("Should the requested output be written to files?")
                    .guiName("Write to files")
                    .build();

    private PluginParameter<String> outputName =
            new PluginParameter.Builder<>("savePath", "", String.class)
                    .description("The base file path for the save files. Each file saved will add a descriptive name to the base name.")
                    .guiName("Base file path")
                    .outFile()
                    .dependentOnParameter(writeFiles)
                    .build();

    public StepwiseAdditiveModelFitterPlugin() {
        super(null, false);
    }

    public StepwiseAdditiveModelFitterPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        //input data should be a single GenotypePhenotype
        List<Datum> datumList = input.getDataOfType(GenotypePhenotype.class);
        if (datumList.size() != 1)
            throw new IllegalArgumentException("Choose exactly one dataset that has combined genotype and phenotype data.");
        myGenoPheno = (GenotypePhenotype) datumList.get(0).getData();
        datasetName = datumList.get(0).getName();

        myFactorNameList =
                myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor).stream()
                        .map(pa -> pa.name())
                        .collect(Collectors.toList());

        myFactorNameList.add(0, "None");
        nestingFactor = PluginParameter.getInstance(nestingFactor, myFactorNameList);
    }

    @Override
    public DataSet processData(DataSet input) {
        StepwiseAdditiveModelFitter stamFitter =
                new StepwiseAdditiveModelFitter(myGenoPheno, datasetName);

        //set parameters
        if (usePermutations.value()) {
            stamFitter.numberOfPermutations(numberOfPermutations.value().intValue());
            stamFitter.permutationAlpha(permutationAlpha.value());
        } else {
            stamFitter.enterLimit(enterLimit.value());
            stamFitter.exitLimit(exitLimit.value());
        }

        if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability)
            stamFitter.useReferenceProbability(true);
        else
            stamFitter.useReferenceProbability(false);

        if (isNested.value()) {
            List nestingList = nestingFactor.value();
            if (nestingFactor.isEmpty()) {
                stamFitter.isNested(false);
            } else {
                String factorName = (String) nestingList.get(0);
                if (factorName.equals("None"))
                    stamFitter.isNested(false);
                else {
                    stamFitter.isNested(true);
                    stamFitter.nestingEffectName(factorName);
                }
            }

        } else {
            stamFitter.isNested(false);
        }

        stamFitter.maxSitesInModel(maxTerms.value());
        stamFitter.modelSelectionCriterion(modelCriterion.value());
        stamFitter.createAnovaReport(createAnova.value());
        stamFitter.createPostScanEffectsReport(createEffects.value());
        stamFitter.createPreScanEffectsReport(createEffectsPrescan.value());
        stamFitter.createResidualsByChr(createResiduals.value());
        stamFitter.createStepReport(createStep.value());

        List<Datum> resultData = new ArrayList<>();
        if (createStep.value()) {
            String name = "";
            String comment = "";
            TableReport output = stamFitter.getSteps();
            resultData.add(new Datum(name, output, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_steps.txt";
                TableReportUtils.saveDelimitedTableReport(output, new File(filename));
            }
        }
        if (createAnova.value()) {
            TableReport output1 = stamFitter.getAnovaReport();
            String name = "";
            String comment = "";
            resultData.add(new Datum(name, output1, comment));

            TableReport output2 = stamFitter.getAnovaReportWithCI();
            name = "";
            comment = "";
            resultData.add(new Datum(name, output2, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_prescan_anova.txt";
                TableReportUtils.saveDelimitedTableReport(output1, new File(filename));
                filename = outputName.value() + "_anova.txt";
                TableReportUtils.saveDelimitedTableReport(output2, new File(filename));
            }
        }
        if (createEffectsPrescan.value()) {
            String name = "";
            String comment = "";
            TableReport output = stamFitter.getMarkerEffectReport();
            resultData.add(new Datum(name, output, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_prescan_effects.txt";
                TableReportUtils.saveDelimitedTableReport(output, new File(filename));
            }
        }
        if (createEffects.value()) {
            String name = "";
            String comment = "";
            TableReport output = stamFitter.getMarkerEffectReportWithCI();
            resultData.add(new Datum(name, output, comment));

            if (writeFiles.value()) {
                String filename = outputName.value() + "_effects.txt";
                TableReportUtils.saveDelimitedTableReport(output, new File(filename));
            }
        }
        if (createResiduals.value()) {
            List<Phenotype> phenoList = stamFitter.getResidualPhenotypesByChromosome();
            phenoList.stream().forEach(pheno -> {
                String name = "";
                String comment = "";
                resultData.add(new Datum(name, pheno, comment));
            });

            if (writeFiles.value()) {
                phenoList.stream().forEach(pheno -> {
                    String traitname = pheno.attributeName(pheno.numberOfAttributes() - 1);
                    String filename = outputName.value() + "_" + traitname + ".txt";
                    PhenotypeUtils.write(pheno, filename);
                });
            }
        }

        return new DataSet(resultData, this);
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = StepwiseOLSModelFitterPlugin.class.getResource("stepwise.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Stepwise-Multithread";
    }

    @Override
    public String getToolTipText() {
        return "Fit a model using stepwise forward-backward regression.";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(StepwiseAdditiveModelFitterPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public DataSet runPlugin(DataSet input) {
        return null;
    }

    /**
     * The model selection criterion used to determine which
     * terms enter the model and how many. Value must be one
     * of pval, bic, mbic, or aic
     *
     * @return Model selection criterion
     */
    public AdditiveSite.CRITERION modelCriterion() {
        return modelCriterion.value();
    }

    /**
     * Set Model selection criterion. The model selection
     * criterion used to determine which terms enter the model
     * and how many. Value must be one of pval, bic, mbic,
     * or aic
     *
     * @param value Model selection criterion
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin modelCriterion(AdditiveSite.CRITERION value) {
        modelCriterion = new PluginParameter<>(modelCriterion, value);
        return this;
    }

    /**
     * Should permutations be used to set the enter and exit
     * limits for stepwise regression? A permutation test
     * will be used to determine the enter limit. The exit
     * limit will be set to 2 times the enter limit.
     *
     * @return Use permutations
     */
    public Boolean usePermutations() {
        return usePermutations.value();
    }

    /**
     * Set Use permutations. Should permutations be used to
     * set the enter and exit limits for stepwise regression?
     * A permutation test will be used to determine the enter
     * limit. The exit limit will be set to 2 times the enter
     * limit.
     *
     * @param value Use permutations
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin usePermutations(Boolean value) {
        usePermutations = new PluginParameter<>(usePermutations, value);
        return this;
    }

    /**
     * The number of permutations used to determine the enter
     * limit.
     *
     * @return Number of permutations
     */
    public Integer numberOfPermutations() {
        return numberOfPermutations.value();
    }

    /**
     * Set Number of permutations. The number of permutations
     * used to determine the enter limit.
     *
     * @param value Number of permutations
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin numberOfPermutations(Integer value) {
        numberOfPermutations = new PluginParameter<>(numberOfPermutations, value);
        return this;
    }

    /**
     * Type I errors will be controlled at this level.
     *
     * @return Alpha for permutations
     */
    public Double permutationAlpha() {
        return permutationAlpha.value();
    }

    /**
     * Set Alpha for permutations. Type I errors will be controlled
     * at this level.
     *
     * @param value Alpha for permutations
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin permutationAlpha(Double value) {
        permutationAlpha = new PluginParameter<>(permutationAlpha, value);
        return this;
    }

    /**
     * When p-value is the model selection criteria, model
     * fitting will stop when the next term chosen has a p-value
     * greater than the enterLimit. This value will be over-ridden
     * the permutation test, if used.
     *
     * @return enterLimit
     */
    public Double enterLimit() {
        return enterLimit.value();
    }

    /**
     * Set enterLimit. When p-value is the model selection
     * criteria, model fitting will stop when the next term
     * chosen has a p-value greater than the enterLimit. This
     * value will be over-ridden the permutation test, if
     * used.
     *
     * @param value enterLimit
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin enterLimit(Double value) {
        enterLimit = new PluginParameter<>(enterLimit, value);
        return this;
    }

    /**
     * During the backward step of model fitting if p-value
     * has been chosen as the selection criterion, if the
     * term in model with the highest p-value has a p-value
     * > exitLimit, it will be removed from the model.
     *
     * @return exitLimit
     */
    public Double exitLimit() {
        return exitLimit.value();
    }

    /**
     * Set exitLimit. During the backward step of model fitting
     * if p-value has been chosen as the selection criterion,
     * if the term in model with the highest p-value has a
     * p-value > exitLimit, it will be removed from the model.
     *
     * @param value exitLimit
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin exitLimit(Double value) {
        exitLimit = new PluginParameter<>(exitLimit, value);
        return this;
    }

    /**
     * The maximum number of SNPs/markers that will be fit.
     * If the model selection criterion is met first, then
     * the fitting process will stop at that point.
     *
     * @return Max SNPs/markers
     */
    public Integer maxTerms() {
        return maxTerms.value();
    }

    /**
     * Set Max SNPs/markers. The maximum number of SNPs/markers
     * that will be fit. If the model selection criterion
     * is met first, then the fitting process will stop at
     * that point.
     *
     * @param value Max SNPs/markers
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin maxTerms(Integer value) {
        maxTerms = new PluginParameter<>(maxTerms, value);
        return this;
    }

    /**
     * Should SNPs/markers be nested within a factor, such
     * as family?
     *
     * @return Is Nested
     */
    public Boolean isNested() {
        return isNested.value();
    }

    /**
     * Set Is Nested. Should SNPs/markers be nested within
     * a factor, such as family?
     *
     * @param value Is Nested
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin isNested(Boolean value) {
        isNested = new PluginParameter<>(isNested, value);
        return this;
    }

    /**
     * Nest markers within this factor.
     *
     * @return Nesting factor
     */
    public List nestingFactor() {
        return nestingFactor.value();
    }

    /**
     * Set Nesting factor. Nest markers within this factor.
     *
     * @param value Nesting factor
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin nestingFactor(List value) {
        nestingFactor = new PluginParameter<>(nestingFactor, value);
        return this;
    }

    /**
     * If the genotype table contains more than one type of
     * genotype data, choose the type to use for the analysis.
     *
     * @return Genotype Component
     */
    public GenotypeTable.GENOTYPE_TABLE_COMPONENT genotypeTable() {
        return myGenotypeTable.value();
    }

    /**
     * Set Genotype Component. If the genotype table contains
     * more than one type of genotype data, choose the type
     * to use for the analysis.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    public StepwiseAdditiveModelFitterPlugin genotypeTable(GenotypeTable.GENOTYPE_TABLE_COMPONENT value) {
        myGenotypeTable = new PluginParameter<>(myGenotypeTable, value);
        return this;
    }

}
