package net.maizegenetics.analysis.association;

import java.util.ArrayList;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastNestedModel;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReportBuilder;

public class DiscreteSitesFullModelOnlyFELM extends AbstractFixedEffectLM {
	protected int numberOfSiteReportColumns;
	protected int numberOfAlleleReportColumns;
	String[] siteGenotypes;
	
	public DiscreteSitesFullModelOnlyFELM(Datum dataset, FixedEffectLMPlugin parentPlugin) {
		super(dataset, parentPlugin);
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
	protected void analyzeSite() {
		
		//add the marker to the model
		myModel = new ArrayList<ModelEffect>(myBaseModel);
        
        ArrayList<String> markerIds = new ArrayList<>();
        int[] markerLevels = ModelEffectUtils.getIntegerLevels(siteGenotypes, markerIds);
        String siteName = myGenoPheno.genotypeTable().siteName(myCurrentSite);
        FactorModelEffect markerEffect = new FactorModelEffect(markerLevels, true, siteName);
        myModel.add(markerEffect);
        
        //add taxa:marker effect at end
        taxaEffectNumber = -1;
        if (areTaxaReplicated) {
        	taxaEffectNumber = myModel.size();
        	myModel.add(taxaEffect());
        }
        
        //calculate model
        double[] modelSSdf;
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, siteData);
        modelSSdf = sflm.getModelcfmSSdf();
        markerSSdf = sflm.getIncrementalSSdf(numberOfBaseEffects);
        if (areTaxaReplicated) errorSSdf = sflm.getIncrementalSSdf(taxaEffectNumber);
        else errorSSdf = sflm.getResidualSSdf();
        double F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
        double p;
        try {
        	p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
        } catch (Exception e) {
        	p = Double.NaN;
        }
        double[] beta = sflm.getBeta();
		G = null;
		if (permute) G = sflm.getInverseOfXtX();
        
        //add results to site report
        //columns = {"Trait","Marker","Chr","Position","marker_F","marker_p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" }
        Object[] rowData = new Object[numberOfSiteReportColumns];
        int columnCount = 0;
        rowData[columnCount++] = currentTraitName;
        rowData[columnCount++] = siteName;
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
        rowData[columnCount++] = new Double(F);
        rowData[columnCount++] = new Double(p);
        if (permute) rowData[columnCount++] = "";
        rowData[columnCount++] = new Double(markerSSdf[0]/(modelSSdf[0] + errorSSdf[0]));
        rowData[columnCount++] = new Double(markerSSdf[1]);
        rowData[columnCount++] = new Double(markerSSdf[0]/markerSSdf[1]);
        rowData[columnCount++] = new Double(errorSSdf[1]);
        rowData[columnCount++] = new Double(errorSSdf[0]/errorSSdf[1]);
        rowData[columnCount++] = new Double(modelSSdf[1]);
        rowData[columnCount++] = new Double(modelSSdf[0]/modelSSdf[1]);
        siteReportBuilder.add(rowData);
        
        //add results to allele report
        //columns = {"Trait","Marker","Chr","Position","Obs","Allele","Estimate"}
        int numberOfAlleles = markerIds.size();
        int[] alleleCounts = markerEffect.getLevelCounts();
        for (int a = 0; a < numberOfAlleles; a++) {
            rowData = new Object[numberOfAlleleReportColumns];
            columnCount = 0;
            rowData[columnCount++] = currentTraitName;
            rowData[columnCount++] = siteName;
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
            rowData[columnCount++] = new Integer(alleleCounts[a]);
            rowData[columnCount++] = markerIds.get(a);
            if (a < numberOfAlleles - 1) rowData[columnCount++] = new Double(beta[numberOfBaseEffects + a]);
            else rowData[columnCount++] = new Double(0.0);
            alleleReportBuilder.add(rowData);
        }
        
	}

	@Override
	protected void getGenotypeAndUpdateMissing(BitSet missingObsBeforeSite) {
		String[] allSiteGenotypes = myGenoPheno.getStringGenotype(myCurrentSite);
		int n = allSiteGenotypes.length;
		missingObsForSite = new OpenBitSet(missingObsBeforeSite);
		for (int i = 0; i < n; i++) {
			if (allSiteGenotypes[i].contains("N")) missingObsForSite.fastSet(i);
		}
		siteGenotypes = AssociationUtils.getNonMissingValues(allSiteGenotypes, missingObsForSite);
	}

}
