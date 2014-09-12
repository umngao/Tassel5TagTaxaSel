package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import java.util.Arrays;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.util.TableReportBuilder;

public class DiscreteSitesAddDomFELM extends AbstractFixedEffectLM {
	protected int numberOfSiteReportColumns;
	protected int numberOfAlleleReportColumns;
	
	public DiscreteSitesAddDomFELM(Datum dataset) {
		super(dataset);
	}

	@Override
	public void initializeReportBuilders() {
		String tableName = "GLM Site Tests - " + myDatum.getName();
		String[] columnNames;
		if (permute) columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","perm_p","marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
		else  columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","marker_Rsq","marker_df","marker_MS","error_df","error_MS","model_df","model_MS" };
		numberOfSiteReportColumns = columnNames.length;
		siteReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		tableName = "GLM Allele Estimates - " + myDatum.getName();
		columnNames = new String[]{"Trait","Marker","Chr","Position","Obs","Allele","Estimate"};
		numberOfAlleleReportColumns = columnNames.length;
		alleleReportBuilder = TableReportBuilder.getInstance(tableName, columnNames);
		siteTableReportRows = new ArrayList<Object[]>();
		markerpvalueColumn = 5;
		permpvalueColumn = 6;
	}

	@Override
	protected void analyzeSite(int siteNumber, ArrayList<ModelEffect> model) {
		myModel = new ArrayList<ModelEffect>(model);
		GenotypeTable myGenotype = myGenoPheno.genotypeTable();
		
		//solve the full model
		//add the marker to the model
        String[] siteGenotypes = getNonMissingValues(myGenoPheno.getStringGenotype(siteNumber), missingObsForSite);
        ArrayList<String> markerIds = new ArrayList<>();
        int[] markerLevels = ModelEffectUtils.getIntegerLevels(siteGenotypes, markerIds);
        String siteName = myGenotype.siteName(siteNumber);
        FactorModelEffect markerEffect = new FactorModelEffect(markerLevels, true, siteName);
        myModel.add(markerEffect);
        
        //add taxa:marker effect at end
        taxaEffectNumber = -1;
        if (areTaxaReplicated) {
        	taxaEffectNumber = myModel.size();
        	myModel.add(taxaEffect());
        }

        //calculate model
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, siteData);
        double[] modelSSdf = sflm.getModelcfmSSdf();
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
 
		//calculate additive and dominance F and p
		ArrayList<ModelEffect> myAdditiveModel = new ArrayList<>(model);
		byte[] alleles = myGenotype.alleles(siteNumber);
		int nAlleles = alleles.length;
		for (int a = 0; a < nAlleles - 1; a++) {
			byte[] allMarker = myGenotype.genotypeAllTaxa(siteNumber);
			byte[] marker = getNonMissingBytes(allMarker, missingObsForSite);
			myAdditiveModel.add(new CovariateModelEffect(ModelEffectUtils.getNumericCodingForAdditiveModel(marker, alleles[a])));
		}
		
        //add taxa:marker effect at end
        int addModelTaxaEffectNumber = -1;
        if (areTaxaReplicated) {
        	addModelTaxaEffectNumber = myAdditiveModel.size();
        	myAdditiveModel.add(taxaEffect());
        }

		
		SweepFastLinearModel sflmAdd = new SweepFastLinearModel(myAdditiveModel, siteData);
		double[] additiveErrorSSdf;
        if (areTaxaReplicated) additiveErrorSSdf = sflmAdd.getIncrementalSSdf(addModelTaxaEffectNumber);
        else additiveErrorSSdf = sflmAdd.getResidualSSdf();

		double[] addTermSSdf = new double[2];
		for (int a = 0; a < nAlleles - 1; a++) {
			double[] thisTermSSdf = sflmAdd.getIncrementalSSdf(a + numberOfBaseEffects);
			addTermSSdf[0] += thisTermSSdf[0];
			addTermSSdf[1] += thisTermSSdf[1];
		}
        double Fadd = addTermSSdf[0] / addTermSSdf[1] / additiveErrorSSdf[0] * additiveErrorSSdf[1];
        double padd;
        try {
            padd = LinearModelUtils.Ftest(Fadd, markerSSdf[1], additiveErrorSSdf[1]);
        } catch (Exception e) {
            padd = Double.NaN;
        }

		//dominance term, F = (reduction in SS Error from full model compared to additive model)/ reduction in df / error MS for full model
        double[] domTermSSdf = new double[2];
        domTermSSdf[0] = markerSSdf[0] - addTermSSdf[0];
        domTermSSdf[1] = markerSSdf[1] - addTermSSdf[1];
        double Fdom = domTermSSdf[0] / domTermSSdf[1] / errorSSdf[0] * errorSSdf[1];
        double pdom;
        try {
        	pdom = LinearModelUtils.Ftest(Fdom, markerSSdf[1], errorSSdf[1]);
        } catch (Exception e) {
        	pdom = Double.NaN;
        }

        //add results to site report
        // column names {"Trait","Marker","Chr","Position","marker_F","marker_p","perm_p","marker_Rsq","add_F","add_p","dom_F","dom_p", "marker_df","marker_MS","error_df","error_MS","model_df","model_MS" }
        Object[] rowData = new Object[numberOfSiteReportColumns];
        int columnCount = 0;
        rowData[columnCount++] = currentTraitName;
        rowData[columnCount++] = siteName;
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(siteNumber);
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(siteNumber);
        rowData[columnCount++] = new Double(F);
        rowData[columnCount++] = new Double(p);
        if (permute) rowData[columnCount++] = "";
        rowData[columnCount++] = new Double(markerSSdf[0]/(modelSSdf[0] + errorSSdf[0]));
        rowData[columnCount++] = new Double(Fadd);
        rowData[columnCount++] = new Double(padd);
        rowData[columnCount++] = new Double(Fdom);
        rowData[columnCount++] = new Double(pdom);
        rowData[columnCount++] = new Double(markerSSdf[1]);
        rowData[columnCount++] = new Double(markerSSdf[0]/markerSSdf[1]);
        rowData[columnCount++] = new Double(errorSSdf[1]);
        rowData[columnCount++] = new Double(errorSSdf[0]/errorSSdf[1]);
        rowData[columnCount++] = new Double(modelSSdf[1]);
        rowData[columnCount++] = new Double(modelSSdf[0]/modelSSdf[1]);
        siteReportBuilder.add(rowData);
        if (permute) siteTableReportRows.add(rowData);
        
        //add results to allele report
        int numberOfAlleles = markerIds.size();
        int[] alleleCounts = markerEffect.getLevelCounts();
        int firstEstimateIndex = beta.length - numberOfAlleles + 1;
        for (int a = 0; a < numberOfAlleles; a++) {
            rowData = new Object[numberOfAlleleReportColumns];
            columnCount = 0;
            rowData[columnCount++] = currentTraitName;
            rowData[columnCount++] = siteName;
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(siteNumber);
            rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(siteNumber);
            rowData[columnCount++] = new Integer(alleleCounts[a]);
            rowData[columnCount++] = markerIds.get(a);
            if (a == numberOfAlleles - 1) rowData[columnCount++] = new Double(beta[firstEstimateIndex + a]);
            else rowData[columnCount++] = new Double(0.0);
            alleleReportBuilder.add(rowData);
        }

	}

}
