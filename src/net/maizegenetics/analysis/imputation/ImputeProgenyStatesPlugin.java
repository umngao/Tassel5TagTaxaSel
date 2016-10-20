package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

public class ImputeProgenyStatesPlugin extends AbstractPlugin {
	
	private PluginParameter<Boolean> rephase = new PluginParameter.Builder<>("rephase", false, Boolean.class)
			.guiName("Rephase Parents First")
			.description("If true, rephase parents before imputing. If false, input haplotypes from file.")
			.build();
	private PluginParameter<String> input = PluginParameter.getLabelInstance("Input --------------------");
	private PluginParameter<String> parentageFile = new PluginParameter.Builder<>("parentage", null, String.class)
			.guiName("Parentage Input File")
			.inFile()
			.description("The input file containing the parentage which lists the parents of each progeny and whether they were derived by self or outcross.")
			.required(true)
			.build();
	private PluginParameter<String> parentHaplotypeFilename = new PluginParameter.Builder<>("parentHap", null, String.class)
			.guiName("Parent Haplotypes Input File")
			.inFile()
			.description("The input file containing the parent haplotypes expressed as nucleotides. Only needed when \"Use Haplotype Probabilities\" (-prob) = false.")
			.build();
	private PluginParameter<String> progenyFile = new PluginParameter.Builder<>("progeny", null, String.class)
			.dependentOnParameter(rephase)
			.guiName("Progeny States Input File")
			.inFile()
			.description("The input file containing the progeny states (parentcalls). Only needed for rephasing using haplotype probabilities.")
			.build();
	private PluginParameter<String> output = PluginParameter.getLabelInstance("Output --------------------");
	private PluginParameter<String> imputedFile = new PluginParameter.Builder<>("imputedOut", null, String.class)
			.guiName("Imputed Genotypes Output File")
			.outFile()
			.description("The output file containing the imputed progeny genotypes in hapmap format.")
			.build();
	private PluginParameter<String> statesFile = new PluginParameter.Builder<>("statesOut", null, String.class)
			.guiName("Progeny States Output File")
			.outFile()
			.description("The output file containing the new progeny states (parentcalls) in hapmap format")
			.build();
	private PluginParameter<String> hapProbFile = new PluginParameter.Builder<>("probOut", null, String.class)
			.guiName("Updated Haplotype Probabilities Output File")
			.outFile()
			.description("The output file containing the new parent haplotype probabilities, binary format. A .bin extension will be appended if not present.")
			.build();

	
	public ImputeProgenyStatesPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	@Override
	public DataSet processData(DataSet input) {
		List<Datum> resultList = new ArrayList<>();
		GenotypeTable inputGenotype = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
		ImputeCrossProgeny icp = new ImputeCrossProgeny();
		icp.setMyGenotype(inputGenotype); //input
		icp.setParentage(parentageFile.value()); //input
		icp.setHaplotypeMap(parentHaplotypeFilename.value());  //input
		icp.setImputedGenotypeOutFilename(imputedFile.value());  //output
		icp.setParentcallOutFilename(statesFile.value());  //output
		icp.setPhasedParentOutFilename(hapProbFile.value());  //output
		
		if (rephase.value()) {
			
			icp.setParentCallInputFilename(progenyFile.value());
			icp.improveImputedProgenyStates();
		} else {
			icp.imputeAll();
		}
		
		return new DataSet(resultList, this);
	}

	@Override
	public ImageIcon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getButtonName() {
		return "Progeny States";
	}

	@Override
	public String getToolTipText() {
		return "Impute progeny states from parent haplotypes.";
	}

	@Override
	public String pluginDescription() {
		return "Impute progeny states using parental haplotypes either represented as nucleotides or as the probability that "
				+ "a haplotype carries the major allele at a site. The plugin also provides a method for estimating haplotype "
				+ "probabilities from progeny states.";
	}

	
}
