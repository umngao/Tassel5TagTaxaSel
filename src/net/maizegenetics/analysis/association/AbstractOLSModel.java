package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SolveByOrthogonalizing;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public abstract class AbstractOLSModel implements FixedEffectLM {
	protected final Datum myDatum;
	protected final GenotypePhenotype myGenoPheno;
	protected final int numberOfObservations;
	protected final int numberOfSites;
	protected final List<PhenotypeAttribute> myDataAttributes;
	protected final List<PhenotypeAttribute> myFactorAttributes;
	protected final List<PhenotypeAttribute> myCovariateAttributes;
	protected TableReportBuilder siteReportBuilder;
	protected TableReportBuilder alleleReportBuilder;
	protected int numberOfSiteReportColumns;
	protected int numberOfAlleleReportColumns;
	protected float[] allData;
	protected int myCurrentSite;
	protected double[] siteData;
	protected OpenBitSet missingObsForSite;
	protected String currentTraitName;
	protected boolean saveToFile = false;
	protected String siteReportFilename;
	protected String alleleReportFilename;
	protected double maxP = 1.0;
	protected FixedEffectLMPlugin myParentPlugin;
	
	//fields used for permutation testing
	protected boolean permute = false;
	protected int numberOfPermutations = 0;
	protected double[] minP = null;
	protected List<DoubleMatrix> permutedData;
	protected double[] baseErrorSSdf;
	protected double[] totalcfmSSdf;
	protected double[] markerSSdf;
	protected double[] errorSSdf;
	protected List<Object[]> siteTableReportRows;
	protected int markerpvalueColumn;
	protected int permpvalueColumn;
	protected ArrayList<ModelEffect> myModel;
	protected DoubleMatrix G;
	protected ArrayList<ModelEffect> myBaseModel;
	protected int numberOfBaseEffects;
	protected int taxaEffectNumber;
	protected int randomSeed;
	protected boolean useRandomSeed = false;
	protected Random rand = null;


    protected static final Map<SiteScore.SITE_SCORE_TYPE, String> typeNameMap;
    static {
    	typeNameMap = new HashMap<SiteScore.SITE_SCORE_TYPE, String>();
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbA, "A");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbC, "C");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbG, "G");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbT, "T");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbGap, "-");
    	typeNameMap.put(SiteScore.SITE_SCORE_TYPE.ProbInsertion, "+");
    }
    
	public AbstractOLSModel(Datum dataset, FixedEffectLMPlugin parentPlugin) {
		myDatum = dataset;
		myParentPlugin = parentPlugin;
		
		myGenoPheno = (GenotypePhenotype) myDatum.getData();
		numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
		numberOfSites = myGenoPheno.genotypeTable().numberOfSites();
		
		myDataAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.data);
		myFactorAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor);
		myCovariateAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.covariate);
		siteTableReportRows = new ArrayList<Object[]>();
		testTaxaReplication();
	}

	@Override
	public void initializeReportBuilders() {
		String tableName = "GLM Statistics - " + myDatum.getName();
		String[] columnNames = siteReportColumnNames();
		numberOfSiteReportColumns = columnNames.length;
		if (saveToFile) siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames, siteReportFilename);
		else siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		 
		tableName = "GLM Genotype Effects - " + myDatum.getName();
		columnNames = alleleReportColumnNames();
		numberOfAlleleReportColumns = columnNames.length;
		if (saveToFile) alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames, alleleReportFilename);
		else alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
	}
	
	@Override
	public void solve() {
		//loop through data attributes
		//	loop through sites
		initializeReportBuilders();
		int numberOfAttributes = myDataAttributes.size();
		int numberOfTestsTotal = numberOfAttributes * numberOfSites;
		int numberOfTestsCalculated = 0;
		int updateInterval = Math.max(1, numberOfTestsTotal / 100);
		
//		long start = System.currentTimeMillis();
		for (PhenotypeAttribute dataAttribute:myDataAttributes) {
			currentTraitName = dataAttribute.name();
			OpenBitSet missingObs = new OpenBitSet(dataAttribute.missing());
			for (PhenotypeAttribute attr:myFactorAttributes) missingObs.or(attr.missing());
			for (PhenotypeAttribute attr:myCovariateAttributes) missingObs.or(attr.missing());
			allData = (float[]) dataAttribute.allValues();
			if (permute) {
				missingObsForSite = missingObs;
				createPermutedData();
			}
			for (int s = 0; s < numberOfSites; s++) {
				//updata missing obs for this site
				myCurrentSite = s;
				getGenotypeAndUpdateMissing(missingObs);
				siteData = AssociationUtils.getNonMissingDoubles(allData, missingObsForSite);
				myBaseModel = baseModel();
				numberOfBaseEffects = myBaseModel.size();
				analyzeSite();
				if (permute) updateMinP(missingObs);
				
				numberOfTestsCalculated++;
				if (numberOfTestsCalculated % updateInterval == 0) {
				        double percentTested = 100.0 * ((double) numberOfTestsCalculated) / numberOfTestsTotal;
				        percentTested = Math.min(percentTested, 100);
					if (myParentPlugin != null) myParentPlugin.updateProgress((int) percentTested);
				}
			}
//			System.out.printf("Sites analyzed in %d ms\n", System.currentTimeMillis() - start);
			if (permute) updateReportsWithPermutationP();
		}
		if (saveToFile) {
			siteReportBuilder.build();
			alleleReportBuilder.build();
		}
	}
	
	@Override
	public TableReport siteReport() {
		saveToFile = true;
		return siteReportBuilder.build();
	}
	
	@Override
	public TableReport alleleReport() {
		saveToFile = true;
		return alleleReportBuilder.build();
	}
	
	@Override
	public List<Datum> datumList() {
		List<Datum> dataList = new ArrayList<Datum>();
		StringBuilder comment = new StringBuilder();
		comment.append("GLM Output\nStatistical Tests for individual variants.\n");
		comment.append("Input data: " + myDatum.getName()).append("\n");
		dataList.add(new Datum("GLM_Stats_" + myDatum.getName(), siteReport(), comment.toString()));
		comment = new StringBuilder();
		comment.append("GLM Output\nGenotype Effect Estimates\n");
		comment.append("Input data: " + myDatum.getName()).append("\n");
		dataList.add(new Datum("GLM_Genotypes_" + myDatum.getName(), alleleReport(), comment.toString()));
		return dataList;
	}
	
	@Override
	public void permutationTest(boolean permute, int nperm) {
		this.permute = permute;
		numberOfPermutations = nperm;
	}

	/**
	 * @param missingObsBeforeSite	a BitSet with bits set for observations missing in model covariates and data
	 */
	protected abstract void getGenotypeAndUpdateMissing(BitSet missingObsBeforeSite);
	
	/**
	 * @param siteNumber		a site number
	 * This method tests the significance of this site and estimates the allele effects then appends the results to the site and allele reports.
	 */
	protected abstract void analyzeSite();
	
	protected String[] siteReportColumnNames() {
		markerpvalueColumn = 5;
		permpvalueColumn = 6;
		if (permute) return new String[]{AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"marker_F",AssociationConstants.STATS_HEADER_P_VALUE,"perm_p","marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
		return new String[] {AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"marker_F",AssociationConstants.STATS_HEADER_P_VALUE,"marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
	}
	
	protected String[] alleleReportColumnNames() {
		return new String[]{AssociationConstants.STATS_HEADER_TRAIT,AssociationConstants.STATS_HEADER_MARKER,AssociationConstants.STATS_HEADER_CHR,AssociationConstants.STATS_HEADER_POSITION,"Obs","Allele","Estimate"};
	}
	
	protected ArrayList<ModelEffect> baseModel() {
		int numberOfNonmissingObs = numberOfObservations - (int) missingObsForSite.cardinality();
		ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
		FactorModelEffect meanEffect = new FactorModelEffect(new int[numberOfNonmissingObs], false);
		meanEffect.setID("mean");
		modelEffects.add(meanEffect);
		
		//add factors to model
		for (PhenotypeAttribute attr:myFactorAttributes) {
			String[] factorLabels = AssociationUtils.getNonMissingValues((String[]) attr.allValues(), missingObsForSite);
			FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, attr.name());
			modelEffects.add(fme);
		}

		//add covariates to model
		for (PhenotypeAttribute attr:myCovariateAttributes) {
			double[] values = AssociationUtils.getNonMissingDoubles((float[]) attr.allValues(), missingObsForSite);
			CovariateModelEffect cme = new CovariateModelEffect(values, attr.name());
			modelEffects.add(cme);
		}

		return modelEffects;
	}
	
	protected void createPermutedData() {
		permutedData = new LinkedList<>();
		double[] y = AssociationUtils.getNonMissingDoubles(allData, missingObsForSite);
		SweepFastLinearModel sflm = new SweepFastLinearModel(baseModel(), y);
		DoubleMatrix residuals = sflm.getResiduals();
		DoubleMatrix predicted = sflm.getPredictedValues();
		baseErrorSSdf = sflm.getResidualSSdf();
		totalcfmSSdf = new double[2];
		double[] modelSSdf = sflm.getModelcfmSSdf();
		totalcfmSSdf[0] = baseErrorSSdf[0] + modelSSdf[0];
		totalcfmSSdf[1] = baseErrorSSdf[1] + modelSSdf[1];
		
		if (rand == null) {
			if (useRandomSeed) rand = new Random(randomSeed);
			else rand = new Random();
		}
		
		for (int p = 0; p < numberOfPermutations; p++) {
			LinearModelUtils.shuffle(residuals, rand);
			DoubleMatrix permdm = predicted.plus(residuals);
			permutedData.add(permdm);
		}
		
		minP = new double[numberOfPermutations];
		Arrays.fill(minP, 1.0);
	}

	protected void updateReportsWithPermutationP() {
		Arrays.sort(minP);
		for (Object[] row : siteTableReportRows) {
			double pval = ((Double) row[markerpvalueColumn]);
			int ndx = Arrays.binarySearch(minP, pval);
			
			if (ndx < 0) ndx = -(ndx + 1);
			double permPval = (double) (ndx + 1) / (double) numberOfPermutations;
			if (permPval > 1) permPval = 1;
			row[permpvalueColumn] = new Double(permPval);
			
		}
	}

	protected ModelEffect taxaEffect() {
        Taxon[] myTaxa = myGenoPheno.phenotype().taxaAttribute().allTaxa();
        Taxon[] myNonMissingTaxa = AssociationUtils.getNonMissingValues(myTaxa, missingObsForSite);
        int[] taxaLevels = ModelEffectUtils.getIntegerLevels(myNonMissingTaxa);
        FactorModelEffect taxaEffect = new FactorModelEffect(taxaLevels, true, "Taxon");
        return taxaEffect;
	}
	
	protected void testTaxaReplication() {
		int numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
		int numberOfTaxa = myGenoPheno.genotypeTable().numberOfTaxa();
		String errmsg = String.format("GLM GWAS no longer supports replicated taxa: number of observations = %d, number of taxa = %d.", numberOfObservations, numberOfTaxa);
		if (numberOfTaxa < numberOfObservations) throw new IllegalArgumentException(errmsg);
	}
	
	protected void updateMinP(BitSet missingObsBeforeSite) {
		int numberOfObsTotal = allData.length;
		int sizeOfPermutedData = permutedData.get(0).numberOfRows();
		BitSet newMissing = new OpenBitSet(sizeOfPermutedData);
		int permutedDataIndex = -1;
		//this loop finds the observations in the permuted data that are missing for this site
		//it has to account for the fact that missing data due to the phenotype and model factors have already been
		//removed from the permuted data
		for (int i = 0; i < numberOfObsTotal; i++) {
			if (!missingObsBeforeSite.fastGet(i)) {
				permutedDataIndex++;
				if (missingObsForSite.fastGet(i)) newMissing.fastSet(permutedDataIndex);
			}
		}

		int iter = 0;
		int numberOfModelEffects = myModel.size();
		for (DoubleMatrix pdata : permutedData) {
			double[] y = AssociationUtils.getNonMissingDoubles(pdata.to1DArray(), newMissing);
			SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
			markerSSdf = sflm.getIncrementalSSdf(numberOfModelEffects - 1);
			errorSSdf = sflm.getResidualSSdf();
			double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
				if (minP[iter] > p) minP[iter] = p;
			} catch (Exception e) {
				//do nothing
			}

			iter++;
		}

	}
	
	@Override
	public void maxP(double maxP) {
		this.maxP = maxP;
	}

	@Override
	public void siteReportFilepath(String savefile) {
		saveToFile = true;
		siteReportFilename = savefile;
		
	}

	@Override
	public void alleleReportFilepath(String savefile) {
		saveToFile = true;
		alleleReportFilename = savefile;
	}

	/**
	 * This method is used mainly for testing in order to generate reproducible permutation results. 
	 * If the seed is not set, the current time is used to initialize the random number generator.
	 * @param seed	the seed used to initialize the random number generator used by permutation
	 */
	public void setRandomSeed(int seed) {
		randomSeed = seed;
		useRandomSeed = true;
	}
	
}
