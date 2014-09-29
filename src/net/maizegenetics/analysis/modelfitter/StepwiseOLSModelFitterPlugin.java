package net.maizegenetics.analysis.modelfitter;

import java.awt.Frame;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import net.maizegenetics.util.TableReport;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class StepwiseOLSModelFitterPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(StepwiseOLSModelFitterPlugin.class);
    private double enterlimit = 1e-5;
    private double exitlimit = 2e-5;
    private int numberOfPermutations = 0;
    private double alpha = 0.05;
    private double[] enterlimits = null;
    private double[] exitlimits = null;
    int maxNumberOfMarkers = 1000;
    boolean nestMarkers = false;
    int nestingFactorIndex;
    private StepwiseOLSModelFitter.MODEL_TYPE modelType = StepwiseOLSModelFitter.MODEL_TYPE.mbic;

    
    public StepwiseOLSModelFitterPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
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

//		MarkerPhenotypeAdapter theAdapter;
//		if (datasets.get(0).getData() instanceof MarkerPhenotype) {
//			MarkerPhenotype mp = (MarkerPhenotype) datasets.get(0).getData();
//			theAdapter = new MarkerPhenotypeAdapter(mp);
//		} else {
//			theAdapter = new MarkerPhenotypeAdapter((Phenotype) datasets.get(0).getData());
//		}
		GenotypePhenotype myGenoPheno = (GenotypePhenotype) datasets.get(0).getData();
		Phenotype myPhenotype = myGenoPheno.phenotype();
		if (isInteractive()) {
			List<PhenotypeAttribute> factorList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);
			int nfactors = factorList.size();
			String[] mainEffects;
			if (nfactors == 0) mainEffects = null;
			else {
				mainEffects = new String[nfactors];
				for (int i = 0; i < nfactors; i++) mainEffects[i] = factorList.get(i).name();
			}
			StepwiseOLSModelFitterDialog myDialog = new StepwiseOLSModelFitterDialog(mainEffects, getParentFrame());

			if (myDialog.wasCancelled()) {
				myDialog.dispose();
				fireProgress(100);
				return null;
			}
			
			double[] limits = myDialog.getEnterLimits();
			if (limits.length == 1) enterlimit = limits[0];
			else if (limits.length > 1) enterlimits = limits;

			limits = myDialog.getExitLimits();
			if (limits.length == 1) exitlimit = limits[0];
			else if (limits.length > 1) exitlimits = limits;

			nestMarkers = myDialog.isNested();
            modelType = myDialog.getModelType();

                        numberOfPermutations = myDialog.getNumberOfPermutations();
                        alpha = myDialog.getAlpha();

			if (nestMarkers) {
				nestingFactorIndex = myPhenotype.attributeIndexForName(myDialog.getNestedEffect());
			}

			maxNumberOfMarkers = myDialog.getMaxNumberOfMarkers();
			myDialog.dispose();
		}

		StepwiseOLSModelFitter modelFitter = new StepwiseOLSModelFitter(myGenoPheno, datasets.get(0).getName());
		modelFitter.setEnterlimit(enterlimit);
		modelFitter.setExitlimit(exitlimit);
		modelFitter.setEnterlimits(enterlimits);
		modelFitter.setExitlimits(exitlimits);
		modelFitter.setMaxNumberOfMarkers(maxNumberOfMarkers);
		modelFitter.setNested(nestMarkers);
		modelFitter.setNestingEffectIndex(nestingFactorIndex);
		modelFitter.setModelType(modelType);
                modelFitter.setNumberOfPermutations(numberOfPermutations);
                modelFitter.setAlpha(alpha);
                
		modelFitter.runAnalysis();
		
		TableReport trResults = modelFitter.getAnovaReport();
		TableReport trEffects = modelFitter.getMarkerEffectReport();
		TableReport trResultsAfterCIScan = modelFitter.getAnovaReportWithCI();
		TableReport trEffectsAfterCIScan = modelFitter.getMarkerEffectReportWithCI();
                TableReport trPermPvalues = modelFitter.getPermutationReport();
                
		LinkedList<Datum> datumList = new LinkedList<Datum>();
		if (trResults != null) datumList.add(new Datum("ANOVA_stepwise_" + datasets.get(0).getName(), trResults, "comments"));
		if (trEffects != null) datumList.add(new Datum("Marker_estimates_" + datasets.get(0).getName(), trEffects, "comments"));
		if (trResultsAfterCIScan != null) datumList.add(new Datum("ANOVA_stepwise_After_CI_Scan" + datasets.get(0).getName(), trResultsAfterCIScan, "comments"));
		if (trEffectsAfterCIScan != null) datumList.add(new Datum("Marker_estimates_After_CI_Scan" + datasets.get(0).getName(), trEffectsAfterCIScan, "comments"));
		if (trPermPvalues != null) datumList.add(new Datum("Permuted_Pvalues" + datasets.get(0).getName(), trPermPvalues, "comments"));
		
		DataSet myResult = new DataSet(datumList, this);
		fireDataSetReturned(myResult);
 
		return myResult;
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
		return "Stepwise";
	}

	@Override
	public String getToolTipText() {
		return "Fit multiple markers in a single model (experimental).";
	}

	public void setEnterlimit(double enterlimit) {
		this.enterlimit = enterlimit;
	}

	public void setExitlimit(double exitlimit) {
		this.exitlimit = exitlimit;
	}

	public void setEnterlimits(double[] enterlimits) {
		this.enterlimits = enterlimits;
	}

	public void setExitlimits(double[] exitlimits) {
		this.exitlimits = exitlimits;
	}

	public void setNumberOfPermutations(int numberOfPermutations) {
		this.numberOfPermutations = numberOfPermutations;
	}
 
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public void setMaxNumberOfMarkers(int maxNumberOfMarkers) {
		this.maxNumberOfMarkers = maxNumberOfMarkers;
	}

	public void setNestMarkers(boolean nestMarkers) {
		this.nestMarkers = nestMarkers;
	}

	public void setNestingFactorIndex(int nestingFactorIndex) {
		this.nestingFactorIndex = nestingFactorIndex;
	}

	public String getCitation(){
		String citation = "Written in 2013 by Peter Bradbury and Alex Lipka";
		return citation;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg; i++) {
			if (args[i].equals("-t") || args[i].equalsIgnoreCase("-modeltype")) {
				//enum MODEL_TYPE {pvalue, bic, mbic, aic} (from StepwiseOLSModelFitter)
				String parm = args[++i].toLowerCase();
				if (parm.startsWith("pval")) modelType = StepwiseOLSModelFitter.MODEL_TYPE.pvalue;
				else if (parm.equals("bic")) modelType = StepwiseOLSModelFitter.MODEL_TYPE.bic;
				else if (parm.equals("mbic")) modelType = StepwiseOLSModelFitter.MODEL_TYPE.mbic;
				else if (parm.equals("aic")) modelType = StepwiseOLSModelFitter.MODEL_TYPE.aic;
				else {
					myLogger.error("Invalid value for modelType. Defaulting to mbic.");
				}
			}
			else if (args[i].equals("-e") || args[i].equalsIgnoreCase("-enter")) {
				enterlimit = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-x") || args[i].equalsIgnoreCase("-exit")) {
				exitlimit = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-es") || args[i].equalsIgnoreCase("-enterlimits")) {
				enterlimits = parseDoubles(args[++i], ",");
			}
			else if (args[i].equals("-xs") || args[i].equalsIgnoreCase("-exitlimits")) {
				exitlimits = parseDoubles(args[++i], ",");
			}
			else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-maxmarkers")) {
				maxNumberOfMarkers = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-n") || args[i].equalsIgnoreCase("-nest")) {
				String parm = args[++i];
				if (parm.toUpperCase().startsWith("T")) nestMarkers = true;
				else nestMarkers = false;
			}
			else if (args[i].equals("-f") || args[i].equalsIgnoreCase("-nestIndex")) {
				nestingFactorIndex = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("?") || args[i].equals("-?")) myLogger.error(getUsage());
		}
	}
	
	public static double[] parseDoubles(String parms, String delimiter) {
		Pattern delim = Pattern.compile(delimiter);
		String[] vals = delim.split(parms);
		int n = vals.length;
		double[] dbls = new double[n];
		for (int i = 0; i < n; i++) dbls[i] = Double.parseDouble(vals[i]);
		return dbls;
	}
	
	public String getUsage() {
		StringBuilder sb = new StringBuilder("The StepwiseOLSModelFitterPlugin can take the following parameters:\n");
		sb.append("-t or -modeltype: The model selection criteria used to determine which terms enter the model and how many. Value must be one of pvalue, bic, mbic, or aic. (default = mbic)");
		sb.append("-e or -enter : The enter limit or maximum p-value for which a term can enter the model (default = 1e-5).\n");
		sb.append("-x or -exit : A term exits the model on a backward step if its p-value is greater than this value (default = 2e-5.\n");
		sb.append("-es or -enterlimits : a comma separated list of enter limits, one for each trait. Overides global enter limit.\n");
		sb.append("-xs or -exitlimits : a comma separated list of exit limits, one for each trait. Overides global exit limit.\n");
		sb.append("-m or -maxmarkers : the maximum number of markers that will be fit, if the enter limit is not reached first.\n");
		sb.append("-n or -nest : true or false, should markers be nested within a model factor? (default = false).\n");
		sb.append("-f or -nestIndex : if there is more then one factor in the model and nest = true, the index of the nesting factor.\n");
		sb.append("? : prints this usage information.\n");
		return sb.toString();
	}

}
