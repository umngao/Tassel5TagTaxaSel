package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.IntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SolveByOrthogonalizing;
import net.maizegenetics.util.TableReportBuilder;

public class EqtlAssociationPlugin extends AbstractPlugin {
	private GenotypeTable.GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GenotypeTable.GENOTYPE_TABLE_COMPONENT[]{
	        GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability};

	private Datum myDatum;
	private GenotypePhenotype myGenoPheno;
	private GenotypeTable myGenotype;
	private Phenotype myPhenotype;
	private TableReportBuilder myReportBuilder;
	private SolveByOrthogonalizing orthogonalSolver;
	private List<String> phenotypeNames;
	
	//plugin parameter definitions
	private PluginParameter<Double> maxp = new PluginParameter.Builder<>("MaxPValue", .001, Double.class)
			.guiName("MaxPValue")
			.description("The maximum p-value that will be output by the analysis.")
			.build();
    private PluginParameter<Boolean> addOnly = new PluginParameter.Builder<>("addOnly", false, Boolean.class)
    		.description("Should an additive only model be fit? If true, an additive model will be fit. If false, an additive + dominance model will be fit. Default = false.")
    		.guiName("Save to file")
    		.build();
	private PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myGenotypeTable = new PluginParameter.Builder<>("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.class)
			.genotypeTable()
	        .range(GENOTYPE_COMP)
	        .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
	        .build();
    private PluginParameter<Boolean> saveAsFile = new PluginParameter.Builder<>("saveToFile", false, Boolean.class)
    		.description("Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirements"
    				+ "and the results will NOT be saved to the data tree. Default = false.")
    		.guiName("Save to file")
    		.build();
    private PluginParameter<String> reportFilename = new PluginParameter.Builder<>("outputFile", null, String.class)
    		.outFile()
    		.dependentOnParameter(saveAsFile)
    		.description("The name of the file to which these results will be saved.")
    		.guiName("Output File")
    		.build();

	
	public EqtlAssociationPlugin(Frame parentFrame, boolean isInteractive) {
		super (parentFrame, isInteractive);
	}

	public DataSet processData(DataSet input) {
		List<Datum> datumList = input.getDataOfType(GenotypePhenotype.class);
		if (datumList.size() != 1) throw new IllegalArgumentException("Exactly one data set with genotypes and phenotypes must be selected.");
		myDatum = datumList.get(0);
		myGenoPheno = (GenotypePhenotype) myDatum.getData();
		myGenotype = myGenoPheno.genotypeTable();
		myPhenotype = myGenoPheno.phenotype();
		testMissingDataInTheBaseModel();
		initializeOutput();
		initializeOrthogonalizer();
		
		final int nsites = myGenotype.numberOfSites();
    	if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) {
    		if (addOnly.value()) {
        		IntStream.range(0, nsites).parallel()
    			.mapToObj(s-> orthogonalSolver.solveForR(myGenotype.positions().get(s), additiveSite(s)))
    			.forEach(this::updateOutputWithPvalues);
    		} else {
        		IntStream.range(0, nsites).parallel()
    			.mapToObj(s-> orthogonalSolver.solveForR(myGenotype.positions().get(s), additiveSite(s)))
    			.forEach(this::updateOutputWithPvalues);
    		}
    	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
    		IntStream.range(0, nsites).parallel()
    			.mapToObj(s-> orthogonalSolver.solveForR(myGenotype.positions().get(s), referenceProbabilitiesForSite(s)))
    			.forEach(this::updateOutputWithPvalues);
    		
    	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability) {
    		throw new UnsupportedOperationException("Eqtl analysis of allele probabilities is not supported.");
    		//TODO implement
    	}

    	if (!saveAsFile.value()) {
    		String name = "EqtlReport_" + myDatum.getComment();
    		String comment = "Rapid Eqtl analysis.";
    		Datum outDatum = new Datum(name, myReportBuilder.build(), comment);
    		return new DataSet(outDatum, this);
    	}
    	
		return null;
	}
	
	private void initializeOutput() {
		//output is a TableReport with p-value; site position information: chr, position, id; trait name
		String[] columnNames = new String[]{"Trait", "Site_id", "Chr", "Position", "p-value" };
		String name = "EqtlReport_" + myDatum.getName();
		if (saveAsFile.value()) myReportBuilder = TableReportBuilder.getInstance(name, columnNames, reportFilename.value());
		else myReportBuilder = TableReportBuilder.getInstance(name, columnNames);
	}
	
	private void updateOutputWithPvalues(SolveByOrthogonalizing.Marker markerResult) {
		double maxpval = maxp.value();
		double[] pvalues = markerResult.vector2();
		int npheno = pvalues.length;
		Position pos = markerResult.position();
		IntStream.range(0, npheno).filter(i -> pvalues[i] < maxpval)
			.forEach(i -> myReportBuilder.add(new Object[]{phenotypeNames.get(i), pos.getSNPID(), pos.getChromosome().getName(), pos.getPosition(), pvalues[i]}));
	}
	
	private void initializeOrthogonalizer() {
		List<PhenotypeAttribute> phenotypeList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
		List<PhenotypeAttribute> covariateList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);
		List<PhenotypeAttribute> factorList = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);
		
		//build the model, no mean necessary because it will not be used
		List<ModelEffect> baseModel = new ArrayList<>();
		for (PhenotypeAttribute pa : factorList) {
			CategoricalAttribute ca = (CategoricalAttribute) pa;
			baseModel.add(new FactorModelEffect(ca.allIntValues(), true, ca.name()));
		}
		for (PhenotypeAttribute pa : covariateList) {
			NumericAttribute na = (NumericAttribute) pa;
			CovariateModelEffect cme = new CovariateModelEffect(AssociationUtils.convertFloatArrayToDouble(na.floatValues()), na.name());
			baseModel.add(cme);
		}
		
		List<double[]> dataList = phenotypeList.stream()
				.map(pa -> (float[]) pa.allValues())
				.map(a -> AssociationUtils.convertFloatArrayToDouble(a))
				.collect(Collectors.toList());
		
		phenotypeNames = phenotypeList.stream().map(PhenotypeAttribute::name).collect(Collectors.toList());
		
		orthogonalSolver = SolveByOrthogonalizing.getInstanceFromModel(baseModel, dataList);
	}
	
	private void testMissingDataInTheBaseModel() {
		for (PhenotypeAttribute attr:myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor)) {
			String msg = "There is missing data in the factor " + attr.name();
			throw new IllegalArgumentException(msg);
		}
		for (PhenotypeAttribute attr:myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate)) {
			String msg = "There is missing data in the covariate " + attr.name();
			throw new IllegalArgumentException(msg);
		}
		for (PhenotypeAttribute attr:myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data)) {
			String msg = "There is missing data in the phenotype " + attr.name();
			throw new IllegalArgumentException(msg);
		}
	}
	
	private double[] referenceProbabilitiesForSite(int site) {
		int ntaxa = myGenotype.numberOfTaxa();
		double[] probs = IntStream.range(0,  ntaxa).mapToDouble(t -> myGenotype.referenceProbability(t, site)).toArray();
		return replaceNansWithMean(probs);
	}
	
	private double[] additiveSite(int site) {
		//codes genotypes as homozygous major = 1, homozygous minor = -1, het = 0
		//set missing values to the mean
		int ntaxa = myGenotype.numberOfTaxa();
		byte major = myGenotype.majorAllele(site);
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		
		double[] code =  IntStream.range(0, ntaxa).mapToDouble(t -> {
			byte geno = myGenotype.genotype(t, site);
			if (geno == NN) return Double.NaN;
			byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
			double val = 0;
			if (alleles[0] == major) val += 0.5;
			if (alleles[1] == major) val += 0.5;
			return val;
		}).toArray();
		
		return replaceNansWithMean(code);
	}
	
	private double[] replaceNansWithMean(double[] array) {
		//replaces the Nans in array with the mean and returns array as a convenience
		int n = array.length;
		double mean = Arrays.stream(array).sum() / n;
		for (int i = 0; i < n; i++) if (Double.isNaN(array[i])) array[i] = mean;
		return array;
	}
	
	private List<double[]> additiveDominanceSite(int site) {
		//the first double[] is the additive site
		//the second double[] equals 1 for heterozygous sites, is 0 otherwise
		int ntaxa = myGenotype.numberOfTaxa();
		List<double[]> result = new ArrayList<>();
		result.add(additiveSite(site));
		double[] dom = new double[ntaxa];
		for (int t = 0; t < ntaxa; t++) if (myGenotype.isHeterozygous(t, site)) dom[t] = 1;
		result.add(dom);
		return result;
	}
	
	//abstract plugin methods that need to be overridden
	@Override
	public ImageIcon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getButtonName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getToolTipText() {
		// TODO Auto-generated method stub
		return null;
	}
	
	
}
