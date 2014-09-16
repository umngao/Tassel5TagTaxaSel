package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public abstract class AbstractFixedEffectLM implements FixedEffectLM {
	protected final Datum myDatum;
	protected final GenotypePhenotype myGenoPheno;
	protected final int numberOfObservations;
	protected final int numberOfSites;
	protected final List<PhenotypeAttribute> myDataAttributes;
	protected final List<PhenotypeAttribute> myFactorAttributes;
	protected final List<PhenotypeAttribute> myCovariateAttributes;
	protected TableReportBuilder siteReportBuilder;
	protected TableReportBuilder alleleReportBuilder;
	protected float[] allData;
	protected int myCurrentSite;
	protected double[] siteData;
	protected OpenBitSet missingObsForSite;
	protected String currentTraitName;
	protected boolean areTaxaReplicated;
	
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
    
	public AbstractFixedEffectLM(Datum dataset) {
		myDatum = dataset;
		myGenoPheno = (GenotypePhenotype) myDatum.getData();
		numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
		numberOfSites = myGenoPheno.genotypeTable().numberOfSites();
		
		myDataAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.data);
		myFactorAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.factor);
		myCovariateAttributes = myGenoPheno.phenotype().attributeListOfType(ATTRIBUTE_TYPE.covariate);
		initializeReportBuilders();
		siteTableReportRows = new ArrayList<Object[]>();
		testTaxaReplication();
	}

	@Override
	public void initializeReportBuilders() {
		String tableName = "GLM Site Tests - " + myDatum.getName();
		String[] columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
		siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		tableName = "GLM Allele Estimates - " + myDatum.getName();
		columnNames = new String[]{"Trait","Marker","Chr","Position","Obs","Allele","Estimate"};
		alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
	}
	
	@Override
	public void solve() {
		//loop through data attributes
		//	loop through sites
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
				siteData = getNonMissingDoubles(allData, missingObsForSite);
				myBaseModel = baseModel();
				numberOfBaseEffects = myBaseModel.size();
				analyzeSite();
				if (permute) updateMinP();
			}
			
			if (permute) updateReportsWithPermutationP();
		}
	}
	
	@Override
	public TableReport siteReport() {
		return siteReportBuilder.build();
	}
	
	@Override
	public TableReport alleleReport() {
		return alleleReportBuilder.build();
	}
	
	@Override
	public List<Datum> datumList() {
		List<Datum> dataList = new ArrayList<Datum>();
		dataList.add(new Datum("", siteReport(), ""));
		dataList.add(new Datum("", alleleReport(), ""));
		return dataList;
	}
	
	@Override
	public void permutationTest(boolean permute) {
		this.permute = permute;
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
	
	protected ArrayList<ModelEffect> baseModel() {
		int numberOfNonmissingObs = numberOfObservations - (int) missingObsForSite.cardinality();
		ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
		FactorModelEffect meanEffect = new FactorModelEffect(new int[numberOfNonmissingObs], false);
		meanEffect.setID("mean");
		modelEffects.add(meanEffect);
		
		//add factors to model
		for (PhenotypeAttribute attr:myFactorAttributes) {
			String[] factorLabels = getNonMissingValues((String[]) attr.allValues(), missingObsForSite);
			FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, attr.name());
			modelEffects.add(fme);
		}

		//add covariates to model
		for (PhenotypeAttribute attr:myCovariateAttributes) {
			double[] values = getNonMissingDoubles((double[]) attr.allValues(), missingObsForSite);
			CovariateModelEffect cme = new CovariateModelEffect(values, attr.name());
			modelEffects.add(cme);
		}

		return modelEffects;
	}
	
	protected void createPermutedData() {
		permutedData = new LinkedList<>();
		double[] y = getNonMissingDoubles(allData, missingObsForSite);
		SweepFastLinearModel sflm = new SweepFastLinearModel(baseModel(), y);
		DoubleMatrix residuals = sflm.getResiduals();
		DoubleMatrix predicted = sflm.getPredictedValues();
		baseErrorSSdf = sflm.getResidualSSdf();
		totalcfmSSdf = new double[2];
		double[] modelSSdf = sflm.getModelcfmSSdf();
		totalcfmSSdf[0] = baseErrorSSdf[0] + modelSSdf[0];
		totalcfmSSdf[1] = baseErrorSSdf[1] + modelSSdf[1];
		
		Random rand = new Random();
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
			row[permpvalueColumn] = new Double((double) (ndx + 1) / numberOfPermutations);
		}
	}

	protected ModelEffect taxaEffect() {
        Taxon[] myTaxa = myGenoPheno.phenotype().taxaAttribute().allTaxa();
        Taxon[] myNonMissingTaxa = getNonMissingValues(myTaxa, missingObsForSite);
        int[] taxaLevels = ModelEffectUtils.getIntegerLevels(myNonMissingTaxa);
        FactorModelEffect taxaEffect = new FactorModelEffect(taxaLevels, true, "Taxon");
        return taxaEffect;
	}
	
	protected void testTaxaReplication() {
		int numberOfObservations = myGenoPheno.phenotype().numberOfObservations();
		int numberOfTaxa = myGenoPheno.genotypeTable().numberOfTaxa();
		if (numberOfTaxa < numberOfObservations) areTaxaReplicated = true;
		else areTaxaReplicated = false;
	}
	
	protected void updateMinP() {
		if (areTaxaReplicated) {
			//if taxa are replicated
			int iter = 0;
			for (DoubleMatrix pdata : permutedData) {
		        double[] modelSSdf;
		        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, siteData);
		        modelSSdf = sflm.getModelcfmSSdf();
		        markerSSdf = sflm.getIncrementalSSdf(numberOfBaseEffects);
		        errorSSdf = sflm.getIncrementalSSdf(taxaEffectNumber);
		        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
		        double p;
		        try {
		        	p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
		        	if (minP[iter] > p) minP[iter] = p;
		        } catch (Exception e) {
		        	//do nothing
		        }
			}
		} else {
			int iter = 0;
			for (DoubleMatrix pdata : permutedData) {
				DoubleMatrix thisData = getNonMissingValues(pdata, missingObsForSite);
				double yty = pdata.crossproduct(thisData).get(0, 0);
				DoubleMatrix Xty = myModel.get(0).getX().crossproduct(thisData);
				double ssmodel = Xty.crossproduct(G).mult(Xty).get(0, 0);
				double sse = yty - ssmodel;
				double ssmarker = baseErrorSSdf[0] - sse;
				double permF = ssmarker / markerSSdf[1] / sse * errorSSdf[1];
				double permP;
				try {
					permP = LinearModelUtils.Ftest(permF, markerSSdf[1], errorSSdf[1]);
					if (minP[iter] > permP) minP[iter] = permP;
				} catch (Exception e) {
					//do nothing
				}
				iter++;
			}
		}
		
	}

	//protected static methods
	protected static double[] getNonMissingDoubles(double[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		double[] result = new double[resultLength];
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}

	protected static double[] getNonMissingDoubles(float[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		double[] result = new double[resultLength];
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}

	protected static byte[] getNonMissingBytes(byte[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		byte[] result = new byte[resultLength];
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}
	
	protected static <T> T[] getNonMissingValues(T[] allData, BitSet missing) {
		int originalLength = allData.length;
		int resultLength = originalLength - (int) missing.cardinality();
		T[] result = Arrays.copyOf(allData, resultLength);
		int resultCount = 0;
		for (int i = 0; i < originalLength;i++) {
			if (!missing.fastGet(i)) result[resultCount++] = allData[i];
		}
		return result;
	}
	
	protected static DoubleMatrix getNonMissingValues(DoubleMatrix allData, BitSet missing) {
		int originalLength = allData.numberOfRows();
		int numberMissing = (int) missing.cardinality();
		if (numberMissing == 0) return allData;
		if (originalLength > 1) {
			int[] select = new int[originalLength - numberMissing];
			int notMissingCount = 0;
			for (int i = 0; i < originalLength; i++) {
				if (!missing.fastGet(i)) select[notMissingCount++] = i;
			}
			return allData.getSelection(select, null);
		} else {
			originalLength = allData.numberOfColumns();
			int[] select = new int[originalLength - numberMissing];
			int notMissingCount = 0;
			for (int i = 0; i < originalLength; i++) {
				if (!missing.fastGet(i)) select[notMissingCount++] = i;
			}
			return allData.getSelection(null, select);
		}
	}
	
	protected static boolean isMonomorphic(float[] probs) {
		float first = Float.NaN;
		int n = probs.length;
		for (int i = 0; i < n; i++) {
			if (!Float.isNaN(probs[i])) {
				if (Float.isNaN(first)) first = probs[i];
				else if (probs[i] != first) return false;
			}
		}
		return true;
	}
	
}
