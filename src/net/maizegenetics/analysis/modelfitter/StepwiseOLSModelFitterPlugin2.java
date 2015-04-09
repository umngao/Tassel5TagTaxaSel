package net.maizegenetics.analysis.modelfitter;

import com.google.common.collect.Range;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.util.TableReport;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

/**
 * Stepwise Ordinary Least Squares model fitter
 *
 * @author Alex Lipka
 * @author Peter Bradbury
 */
public class StepwiseOLSModelFitterPlugin2 extends AbstractPlugin {

    private PluginParameter<StepwiseOLSModelFitter.MODEL_TYPE> modelType = new PluginParameter.Builder<>("t", StepwiseOLSModelFitter.MODEL_TYPE.bic, StepwiseOLSModelFitter.MODEL_TYPE.class)
            .range(StepwiseOLSModelFitter.MODEL_TYPE.values())
            .guiName("Model type")
            .description("The model selection criteria used to determine which terms enter the model and how many. Value must be one of pvalue, bic, mbic, or aic")
            .build();
    private PluginParameter<Double> enterlimit = new PluginParameter.Builder<>("e", 1e-5, Double.class)
            .range(Range.closed(0.0, 1.0))
            .guiName("Entry limit")
            .description("The enter limit or maximum p-value for which a term can enter the model")
            .build();
    private PluginParameter<Double> exitlimit = new PluginParameter.Builder<>("x", 1e-6, Double.class)
            .range(Range.closed(0.0, 1.0))
            .guiName("Exit limit")
            .description("A term exits the model on a backward step if its p-value is greater than this value")
            .build();
    private PluginParameter<Integer> maxNumberOfMarkers = new PluginParameter.Builder<>("m", 100, Integer.class)
            .range(Range.closed(0, 10000))
            .guiName("Maximum markers")
            .description("The maximum number of markers that will be fit, if the enter limit is not reached first")
            .build();
    private PluginParameter<Boolean> nestMarkers = new PluginParameter.Builder<>("n", false, Boolean.class)
            .guiName("Nest markers")
            .description("Should markers be nested within a model factor")
            .build();
    private PluginParameter<Integer> nestingFactorIndex = new PluginParameter.Builder<>("m", 0, Integer.class)
            .guiName("Nesting factor index")
            .description("If there is more then one factor in the model and nest = true, the index of the nesting factor")
            .build();
    private PluginParameter<Integer> numberOfPermutations = new PluginParameter.Builder<>("m", 0, Integer.class)
            .range(Range.closed(0, 100000))
            .guiName("Number of permutations")
            .description("Number of permutations for the model to determine an empirical alpha")
            .build();


    private static final Logger myLogger = Logger.getLogger(StepwiseOLSModelFitterPlugin2.class);
    private double alpha = 0.05;
    //TODO need to change this to a list
   // private StepwiseOLSModelFitter.MODEL_TYPE modelType = StepwiseOLSModelFitter.MODEL_TYPE.mbic;


    public StepwiseOLSModelFitterPlugin2(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public StepwiseOLSModelFitterPlugin2() {
        super(null, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        List<Datum> datasets = input.getDataOfType(new Class[]{GenotypePhenotype.class});
        if (datasets.size() < 1) {
            String msg = "Error in performFunction: No appropriate dataset selected.";
            myLogger.error(msg);
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), msg, "Error in Model Fitter", JOptionPane.ERROR_MESSAGE);
            }
            return null;
        }

        //only analyze the first dataset
        if (datasets.size() > 1) {
            String msg = "Multiple datasets selected. Only the first will be analyzed.";
            myLogger.info(msg);
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), msg, "Error in Model Fitter", JOptionPane.INFORMATION_MESSAGE);
            }
            return null;
        }

        GenotypePhenotype myGenoPheno = (GenotypePhenotype) datasets.get(0).getData();
        if (nestMarkers.value()) {
            boolean factorCorrect = myGenoPheno.phenotype().attribute(nestingFactorIndex.value()).isTypeCompatible(ATTRIBUTE_TYPE.factor);
            if (!factorCorrect) {
                myLogger.error("Error: Selected factor index for nesting is not a valid factor! Please choose from the list below");
                int[] factorIndices = myGenoPheno.phenotype().attributeIndicesOfType(ATTRIBUTE_TYPE.factor);
                Phenotype mypheno = myGenoPheno.phenotype();
                for (int i : factorIndices) {
                    myLogger.info("     Factor \"" + mypheno.attributeName(i) + "\" is at index " + i);
                }
            }
        }

        StepwiseOLSModelFitter modelFitter = new StepwiseOLSModelFitter(myGenoPheno, datasets.get(0).getName());
        modelFitter.setEnterlimit(enterlimit.value());
        modelFitter.setExitlimit(exitlimit.value());
        modelFitter.setMaxNumberOfMarkers(maxNumberOfMarkers.value());
        modelFitter.setNested(nestMarkers.value());
        modelFitter.setNestingEffectIndex(nestingFactorIndex.value());
        modelFitter.setModelType(modelType.value());
        modelFitter.setNumberOfPermutations(numberOfPermutations.value());
        modelFitter.setAlpha(alpha);

        modelFitter.runAnalysis();

        TableReport trResults = modelFitter.getAnovaReport();
        TableReport trEffects = modelFitter.getMarkerEffectReport();
        TableReport trResultsAfterCIScan = modelFitter.getAnovaReportWithCI();
        TableReport trEffectsAfterCIScan = modelFitter.getMarkerEffectReportWithCI();
        TableReport trPermPvalues = modelFitter.getPermutationReport();

        LinkedList<Datum> datumList = new LinkedList<Datum>();
        if (trResults != null)
            datumList.add(new Datum("ANOVA_stepwise_" + datasets.get(0).getName(), trResults, "comments"));
        if (trEffects != null)
            datumList.add(new Datum("Marker_estimates_" + datasets.get(0).getName(), trEffects, "comments"));
        if (trResultsAfterCIScan != null)
            datumList.add(new Datum("ANOVA_stepwise_After_CI_Scan" + datasets.get(0).getName(), trResultsAfterCIScan, "comments"));
        if (trEffectsAfterCIScan != null)
            datumList.add(new Datum("Marker_estimates_After_CI_Scan" + datasets.get(0).getName(), trEffectsAfterCIScan, "comments"));
        if (trPermPvalues != null)
            datumList.add(new Datum("Permuted_Pvalues" + datasets.get(0).getName(), trPermPvalues, "comments"));

        DataSet myResult = new DataSet(datumList, this);
        fireDataSetReturned(myResult);

        return myResult;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = StepwiseOLSModelFitterPlugin2.class.getResource("stepwise.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Stepwise";
    }

    @Override
    public String getToolTipText() {
        return "Fit multiple markers in a single model (experimental).";
    }

    public String getCitation() {
        String citation = "Written in 2013 by Peter Bradbury and Alex Lipka";
        return citation;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    public static void main(String[] args) {
         GeneratePluginCode.generate(StepwiseOLSModelFitterPlugin2.class);
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

    /**
     * The model selection criteria used to determine which
     * terms enter the model and how many. Value must be one
     * of pvalue, bic, mbic, or aic
     *
     * @return Model type
     */
    public StepwiseOLSModelFitter.MODEL_TYPE modelType() {
        return modelType.value();
    }

    /**
     * Set Model type. The model selection criteria used to
     * determine which terms enter the model and how many.
     * Value must be one of pvalue, bic, mbic, or aic
     *
     * @param value Model type
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 modelType(StepwiseOLSModelFitter.MODEL_TYPE value) {
        modelType = new PluginParameter<>(modelType, value);
        return this;
    }

    /**
     * The enter limit or maximum p-value for which a term
     * can enter the model
     *
     * @return Entry limit
     */
    public Double enterlimit() {
        return enterlimit.value();
    }

    /**
     * Set Entry limit. The enter limit or maximum p-value
     * for which a term can enter the model
     *
     * @param value Entry limit
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 enterlimit(Double value) {
        enterlimit = new PluginParameter<>(enterlimit, value);
        return this;
    }

    /**
     * A term exits the model on a backward step if its p-value
     * is greater than this value
     *
     * @return Exit limit
     */
    public Double exitlimit() {
        return exitlimit.value();
    }

    /**
     * Set Exit limit. A term exits the model on a backward
     * step if its p-value is greater than this value
     *
     * @param value Exit limit
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 exitlimit(Double value) {
        exitlimit = new PluginParameter<>(exitlimit, value);
        return this;
    }

    /**
     * The maximum number of markers that will be fit, if
     * the enter limit is not reached first
     *
     * @return Maximum markers
     */
    public Integer maxNumberOfMarkers() {
        return maxNumberOfMarkers.value();
    }

    /**
     * Set Maximum markers. The maximum number of markers
     * that will be fit, if the enter limit is not reached
     * first
     *
     * @param value Maximum markers
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 maxNumberOfMarkers(Integer value) {
        maxNumberOfMarkers = new PluginParameter<>(maxNumberOfMarkers, value);
        return this;
    }

    /**
     * Should markers be nested within a model factor
     *
     * @return Nest markers
     */
    public Boolean nestMarkers() {
        return nestMarkers.value();
    }

    /**
     * Set Nest markers. Should markers be nested within a
     * model factor
     *
     * @param value Nest markers
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 nestMarkers(Boolean value) {
        nestMarkers = new PluginParameter<>(nestMarkers, value);
        return this;
    }

    /**
     * If there is more then one factor in the model and nest
     * = true, the index of the nesting factor
     *
     * @return Nesting factor index
     */
    public Integer nestingFactorIndex() {
        return nestingFactorIndex.value();
    }

    /**
     * Set Nesting factor index. If there is more then one
     * factor in the model and nest = true, the index of the
     * nesting factor
     *
     * @param value Nesting factor index
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 nestingFactorIndex(Integer value) {
        nestingFactorIndex = new PluginParameter<>(nestingFactorIndex, value);
        return this;
    }

    /**
     * Number of permutations for the model to determine an
     * empirical alpha
     *
     * @return Number of permutations
     */
    public Integer numberOfPermutations() {
        return numberOfPermutations.value();
    }

    /**
     * Set Number of permutations. Number of permutations
     * for the model to determine an empirical alpha
     *
     * @param value Number of permutations
     *
     * @return this plugin
     */
    public StepwiseOLSModelFitterPlugin2 numberOfPermutations(Integer value) {
        numberOfPermutations = new PluginParameter<>(numberOfPermutations, value);
        return this;
    }

}
