package net.maizegenetics.analysis.association;

import java.util.ArrayList;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastNestedModel;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReportBuilder;

public class DiscreteSitesAddDomWithReplicationFELM extends AbstractFixedEffectLM {
	protected int numberOfSiteReportColumns;
	protected int numberOfAlleleReportColumns;
	protected double[] markerSSdf;
	protected double[] errorSSdf;
	ArrayList<ModelEffect> myModel;
	DoubleMatrix G;
	
	public DiscreteSitesAddDomWithReplicationFELM(Datum dataset) {
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
        Taxon[] myTaxa = myGenoPheno.phenotype().taxaAttribute().allTaxa();
        Taxon[] myNonMissingTaxa = getNonMissingValues(myTaxa, missingObsForSite);
        int[] taxaLevels = ModelEffectUtils.getIntegerLevels(myNonMissingTaxa);
        FactorModelEffect taxaEffect = new FactorModelEffect(taxaLevels, true, "Taxon");
        model.add(taxaEffect);
        
        //calculate model
        SweepFastNestedModel sfnm = new SweepFastNestedModel(myModel, siteData);
        double[] taxaSSdf = sfnm.getTaxaInMarkerSSdf();
        double[] markerSSdf = sfnm.getMarkerSSdf();
        double[] errorSSdf = sfnm.getErrorSSdf();
        double[] modelSSdf = sfnm.getModelcfmSSdf();
        double F = markerSSdf[0] / markerSSdf[1] / taxaSSdf[0] * taxaSSdf[1];
        double p;
        try {
            p = LinearModelUtils.Ftest(F, markerSSdf[1], taxaSSdf[1]);
        } catch (Exception e) {
            p = Double.NaN;
        }
        double[] beta = sfnm.getBeta();
		G = null;
		if (permute) G = sfnm.getInverseOfXtX();
 
		//calculate additive and dominance F and p
		ArrayList<ModelEffect> myAdditiveModel = new ArrayList<>(model);
		int firstMarkerEffect = myAdditiveModel.size();
		byte[] alleles = myGenotype.alleles(siteNumber);
		int nAlleles = alleles.length;
		for (int a = 0; a < nAlleles - 1; a++) {
			byte[] allMarker = myGenotype.genotypeAllTaxa(siteNumber);
			byte[] marker = getNonMissingBytes(allMarker, missingObsForSite);
			myAdditiveModel.add(new CovariateModelEffect(ModelEffectUtils.getNumericCodingForAdditiveModel(marker, alleles[a])));
		}
		int taxaEffectNumber = myAdditiveModel.size();
		model.add(taxaEffect);
		SweepFastLinearModel sflmAdd = new SweepFastLinearModel(myAdditiveModel, siteData);
		
		//additive term, F == MS from additive only model / MS Error from additive only model
		double[] addModelSSdf = sflmAdd.getFullModelSSdf();
		double[] addTermSSdf = new double[]{0,0};
		for (int a = 0; a < nAlleles - 1; a++) {
			double[] incrSSdf = sflmAdd.getIncrementalSSdf(firstMarkerEffect + a);
			addTermSSdf[0] += incrSSdf[0];
			addTermSSdf[1] += incrSSdf[1];
		}
		double[] taxaTermSSdf = sflmAdd.getIncrementalSSdf(taxaEffectNumber);
		
        double Fadd = addTermSSdf[0] / addTermSSdf[1] / taxaTermSSdf[0] * taxaTermSSdf[1];
        double padd;
        try {
            padd = LinearModelUtils.Ftest(Fadd, addTermSSdf[1], taxaTermSSdf[1]);
        } catch (Exception e) {
            padd = Double.NaN;
        }

		//dominance term, F = (reduction in SS Error from full model compared to additive model)/ reduction in df / error MS for full model
        double[] domTermSSdf = new double[2];
        domTermSSdf[0] = modelSSdf[0] - addModelSSdf[0];
        domTermSSdf[1] = modelSSdf[1] - addModelSSdf[1];
        double Fdom = domTermSSdf[0] / domTermSSdf[1] / taxaSSdf[0] * taxaSSdf[1];
        double pdom;
        try {
        	pdom = LinearModelUtils.Ftest(Fdom, domTermSSdf[1], taxaSSdf[1]);
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

	@Override
	protected void updateMinP() {
    	int iter = 0;
    	for (DoubleMatrix pdata : permutedData) {
    		double yty = pdata.crossproduct(pdata).get(0, 0);
    		DoubleMatrix Xty = myModel.get(0).getX().crossproduct(pdata);
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
