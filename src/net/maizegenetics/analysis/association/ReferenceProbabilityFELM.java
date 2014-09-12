package net.maizegenetics.analysis.association;

import java.util.ArrayList;

import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.util.TableReportBuilder;

public class ReferenceProbabilityFELM extends AbstractFixedEffectLM {
	int numberOfSiteReportColumns;
	int numberOfAlleleReportColumns;

	public ReferenceProbabilityFELM(Datum data) {
		super(data);
	}
	
	@Override
	public void initializeReportBuilders() {
		String tableName = "GLM Site Tests - " + myDatum.getName();
		String[] columnNames;
		if (permute) columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","perm_p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
		else columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
		numberOfSiteReportColumns = columnNames.length;
		siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		tableName = "GLM Allele Estimates - " + myDatum.getName();
		columnNames = new String[]{"Trait","Marker","Chr","Position","Obs","Allele","Estimate"};
		numberOfAlleleReportColumns = columnNames.length;
		alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		permpvalueColumn = 6;
		markerpvalueColumn = 5;
	}
	
	@Override
	protected void analyzeSite(int siteNumber, ArrayList<ModelEffect> baseEffects) {
		myModel = new ArrayList<ModelEffect>(baseEffects);
		String siteName = myGenoPheno.genotypeTable().siteName(siteNumber);

		float[] probs = alleleProbsOfType(SITE_SCORE_TYPE.ReferenceProbablity, siteNumber);
		double[] covar = getNonMissingDoubles(probs, missingObsForSite);
		myModel.add(new CovariateModelEffect(covar));
		
		if (areTaxaReplicated) myModel.add(taxaEffect());
		
		//solve the model
		SweepFastLinearModel markerModel = new SweepFastLinearModel(myModel, siteData);

        //calculate model
        double[] modelSSdf = markerModel.getModelcfmSSdf();
        markerSSdf = markerModel.getIncrementalSSdf(numberOfBaseEffects);
        if (areTaxaReplicated) errorSSdf = markerModel.getIncrementalSSdf(numberOfBaseEffects + 1);
        else errorSSdf = markerModel.getResidualSSdf();
        
        double rsq = markerSSdf[0] / (modelSSdf[0] + errorSSdf[0]);
        
        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
        double p;
        try {
            p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
        } catch (Exception e) {
            p = Double.NaN;
        }
        double[] beta = markerModel.getBeta();
		if (permute) G = markerModel.getInverseOfXtX();
		
        //add results to site report
        //{"Trait","Marker","Chr","Position","marker_F","marker_p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" }
        Object[] rowData = new Object[numberOfSiteReportColumns];
        int columnCount = 0;
        rowData[columnCount++] = currentTraitName;
        rowData[columnCount++] = siteName;	
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(siteNumber);
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(siteNumber);
        rowData[columnCount++] = new Double(F);
        rowData[columnCount++] = new Double(p);
        if (permute) rowData[columnCount++] = "";
        rowData[columnCount++] = new Double(rsq);
        rowData[columnCount++] = new Double(markerSSdf[1]);
        rowData[columnCount++] = new Double(markerSSdf[0] / markerSSdf[1]);
        rowData[columnCount++] = new Double(errorSSdf[1]);
        rowData[columnCount++] = new Double(errorSSdf[0] / errorSSdf[1]);
        rowData[columnCount++] = new Double(modelSSdf[1]);
        rowData[columnCount++] = new Double(modelSSdf[0] / modelSSdf[1]);
        
        //add results to allele report
        //{"Trait","Marker","Chr","Position","Estimate"}
        int estimateIndex = beta.length - 1;
        rowData = new Object[numberOfAlleleReportColumns];
        columnCount = 0;
        rowData[columnCount++] = currentTraitName;
        rowData[columnCount++] = siteName;
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(siteNumber);
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(siteNumber);
        rowData[columnCount++] = beta[estimateIndex];
        alleleReportBuilder.add(rowData);
	}

}
