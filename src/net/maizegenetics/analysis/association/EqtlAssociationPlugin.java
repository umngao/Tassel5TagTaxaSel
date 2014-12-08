package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ReferenceProbabilitySpliterator;
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
import net.maizegenetics.util.OpenBitSet;
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
	PluginParameter<Double> maxp = new PluginParameter.Builder<>("MaxPValue", .001, Double.class)
			.guiName("MaxPValue")
			.description("The maximum p-value that will be output by the analysis.")
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
		GenotypePhenotype myGenoPheno = (GenotypePhenotype) myDatum.getData();
		Phenotype myPhenotype = myGenoPheno.phenotype();
		GenotypeTable myGenotype = myGenoPheno.genotypeTable();
		testMissingDataInTheBaseModel();
		initializeOutput();
		initializeOrthogonalizer();
		
    	if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype) {
    		throw new UnsupportedOperationException("Eqtl analysis of genotypes is not supported.");
    		//TODO implement
    	} else if (myGenotypeTable.value() == GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
    		ReferenceProbabilitySpliterator splitter = new ReferenceProbabilitySpliterator(myGenotype, 0, myGenotype.numberOfSites());
    		Stream<ReferenceProbabilitySpliterator.ReferenceProbabilityBySite> refProbStream = StreamSupport.stream(splitter, true);
    		refProbStream.map(p -> orthogonalSolver.solveForR(p.myPosition, p.myValues)).forEach(m -> updateOutputWithPvalues(m));
    		
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
