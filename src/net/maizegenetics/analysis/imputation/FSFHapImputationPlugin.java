package net.maizegenetics.analysis.imputation;

import java.awt.Frame;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

public class FSFHapImputationPlugin extends AbstractPlugin {
	
	//parameters for CallParentAllelesPlugin
	private PluginParameter<String> pedigreeFilename = new PluginParameter.Builder<>("pedigrees", null, String.class)
			.description("the pedigree file name")
			.inFile().required(true).build();
	private PluginParameter<Boolean> useClusterAlgorithm = new PluginParameter.Builder<>("cluster", false, Boolean.class)
			.description("use the cluster algorithm").build();
	private PluginParameter<Boolean> useWindowLD = new PluginParameter.Builder<>("windowLD", false, Boolean.class)
			.description("use the windowLD algorithm").build();
	private PluginParameter<Boolean> useBCFilter = new PluginParameter.Builder<>("bc", true, Boolean.class)
			.description("use the single backcross algorithm").build();
	private PluginParameter<Boolean> useMultipleBCFilter = new PluginParameter.Builder<>("multbc", false, Boolean.class)
			.description("use the multiple backcross algorithm").build();
	private PluginParameter<Double> minMinorAlleleFreq = new PluginParameter.Builder<>("minMaf", 0.1, Double.class)
			.description("filter out sites with less than minimumMinorAlleleFrequency").build();
	private PluginParameter<Integer> windowSize = new PluginParameter.Builder<>("window", 50, Integer.class)
			.description("filter out sites with less than minimumMinorAlleleFrequency").build();
	private PluginParameter<Double> minRforSnps = new PluginParameter.Builder<>("minR", 0.2, Double.class)
			.description("filter out sites not correlated with neighboring sites").build();
	private PluginParameter<Double> maxMissing = new PluginParameter.Builder<>("maxMissing", 0.2, Double.class)
			.description("filter out sites with proportion missing > maxMissing").build();
//	private PluginParameter<Boolean> checkSubPops = new PluginParameter.Builder<>("subpops", false, Boolean.class)
//			.description("check subpopulations (rarely used)").build();
	private PluginParameter<Boolean> noHets = new PluginParameter.Builder<>("nohets", false, Boolean.class)
			.description("delete heterozygous calls before imputing").build();
	private PluginParameter<Integer> maxDifference = new PluginParameter.Builder<>("maxDiff", 0, Integer.class)
			.description("use to decide if two haplotypes are equivalent").build();
	private PluginParameter<Integer> minHaplotypeCluster = new PluginParameter.Builder<>("minHap", 5, Integer.class)
			.description("haplotype must be observed at least this often").build();
//	private PluginParameter<Double> maxHetDev = new PluginParameter.Builder<>("maxHetDev", 25.0, Double.class)
//			.description("sites with percent het this many sd above mean not used").build();
	private PluginParameter<Integer> overlap = new PluginParameter.Builder<>("overlap", 25, Integer.class)
			.description("overlap between adjacent windows").build();
	
	//parameters for ViterbiAlgorithmPlugin
	
	//parameters for WritePopulationAlignmentPlugin
	
	public FSFHapImputationPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	@Override
	public String getCitation() {
        return "Swarts K, Li H, Romero Navarro JA, Romay-Alvarez MC, Hearne S, Acharya C, "
                + "Glaubitz JC, Mitchell S, Elshire RJ, Buckler ES, Bradbury PJ (2014) "
                + "FSFHap (Full-Sib Family Haplotype Imputation) and FILLIN "
                + "(Fast, Inbred Line Library ImputatioN) optimize genotypic imputation "
                + "for low-coverage, next-generation sequence data in crop plants. "
                + "Plant Genome (in review)";
	}
	
	@Override
	public String pluginDescription() {
		// TODO Auto-generated method stub
		return "The FSFHapImputation Plugin infers parental haplotypes for a full sib family then uses those haplotypes in an HMM to impute variants. "
				+ "It is effective at correctly imputing heterzygotes in GBS data.";
	}

	public DataSet processData(DataSet input) {
		return null;
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		// TODO Auto-generated method stub
		return "Impute By FSFHap";
	}

	@Override
	public String getToolTipText() {
		return "Impute variants in full sib families";
	}

}
