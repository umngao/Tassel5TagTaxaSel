package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.net.URL;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.log4j.Logger;

import com.google.common.collect.Range;


public class FixedEffectLMPlugin extends AbstractPlugin {
	
    private static final Logger myLogger = Logger.getLogger(FixedEffectLMPlugin.class);
    
	enum GENOTYPE_DATA_TYPE { genotype, probability, allele_probabilities, none };
	private GenotypeTable.GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GenotypeTable.GENOTYPE_TABLE_COMPONENT[]{
	        GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability};
	
    //parameters
	private PluginParameter<Boolean> phenotypeOnly = new PluginParameter.Builder<>("phenoOnly", false, Boolean.class)
			.description("Should the phenotype be analyzed with no markers and BLUEs generated? (BLUE = best linear unbiased estimate)")
			.guiName("Analyze Phenotype Only")
			.build();
    private PluginParameter<Boolean> saveToFile = new PluginParameter.Builder<>("saveToFile", false, Boolean.class)
    		.description("Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirements"
    				+ "and the results will NOT be saved to the data tree. Default = false.")
    		.guiName("Save to file?")
    		.build();
    private PluginParameter<String> siteReportFilename = new PluginParameter.Builder<>("siteFile", null, String.class)
    		.outFile()
    		.dependentOnParameter(saveToFile)
    		.description("The name of the file to which these results will be saved.")
    		.guiName("Save filename")
    		.build();
    private PluginParameter<String> alleleReportFilename = new PluginParameter.Builder<>("alleleFile", null, String.class)
    		.outFile()
    		.dependentOnParameter(saveToFile)
    		.description("The name of the file to which these results will be saved.")
    		.guiName("Save filename")
    		.build();
    private PluginParameter<Double> maxPvalue = new PluginParameter.Builder<>("maxP", 1.0, Double.class)
    		.description("Only results with p <= maxPvalue will be reported. Default = 1.0.")
    		.dependentOnParameter(phenotypeOnly, false)
    		.range(Range.closed(0.0, 1.0))
    		.guiName("max P value")
    		.build();
    private PluginParameter<Boolean> permute = new PluginParameter.Builder<>("permute", false, Boolean.class)
    		.description("Should a permutation analysis be run? The permutation analysis controls the experiment-wise error rate for individual phenotypes.")
    		.dependentOnParameter(phenotypeOnly, false)
    		.guiName("Run Permutations")
    		.build();
    private PluginParameter<Integer> numberOfPermutations = new PluginParameter.Builder<>("nperm", 0, Integer.class)
    		.description("The number of permutations to be run for the permutation analysis.")
    		.dependentOnParameter(permute)
    		.guiName("Number of Permutations")
    		.build();
	private PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myGenotypeTable = new PluginParameter.Builder<>("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.class)
			.genotypeTable()
	        .range(GENOTYPE_COMP)
	        .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
	        .build();
    
    public FixedEffectLMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public String getButtonName() {
        return "GLM";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = FixedEffectLMPlugin.class.getResource("/net/maizegenetics/analysis/images/LinearAssociation.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getToolTipText() {
        return "Use fixed effect model to test associations";
    }

    protected void preProcessParameters(DataSet data) {
    	List<Datum> genoPhenoList = data.getDataOfType(GenotypePhenotype.class);
    	
    	if (genoPhenoList.size() == 0){
    		List<Datum> phenoList = data.getDataOfType(Phenotype.class);
    		if (phenoList.size() == 0) throw new IllegalArgumentException("A dataset that can be analyzed by GLM has not been selected.");
    		else if (phenoList.size() == 1) {
        		phenotypeOnly = new PluginParameter.Builder<>("phenoOnly", true, Boolean.class)
            			.description("Should the phenotype be analyzed with no markers and BLUEs generated? (BLUE = best linear unbiased estimate)")
            			.guiName("Analyze Phenotype Only")
            			.build();
    		} else  throw new IllegalArgumentException("GLM can only process one data set at a time.");
    	} 
    	else if (genoPhenoList.size() > 1)  throw new IllegalArgumentException("GLM can only process one data set at a time.");
    			
    }
    
    public DataSet processData(DataSet data) {
    	if (phenotypeOnly.value()) {
    		Datum myDatum = data.getDataOfType(Phenotype.class).get(0);
    		PhenotypeLM plm = new PhenotypeLM(myDatum);
    		return new DataSet(plm.datumList(), this);
    	} else {
        	FixedEffectLM myLM; 
    		Datum myDatum = data.getDataOfType(GenotypePhenotype.class).get(0);
        	if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) {
        		myLM = new DiscreteSitesFELM(myDatum);
        	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
        		myLM = new ReferenceProbabilityFELM(myDatum);
        	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability) {
        		myLM = new AlleleProbabilityFELM(myDatum);
        	} else return null;
        	if (permute.value()) myLM.permutationTest(true, numberOfPermutations.value());
        	if (saveToFile.value()) {
        		myLM.siteReportFilepath(siteReportFilename.value());
        		myLM.alleleReportFilepath(alleleReportFilename.value());
        	}
        	myLM.maxP(maxPvalue.value());
        	myLM.solve();
        	if (saveToFile.value()) return null;
        	else return new DataSet(myLM.datumList(), this);
    	} 
    	
    }
}
