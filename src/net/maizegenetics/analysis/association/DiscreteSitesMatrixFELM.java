package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import java.util.Arrays;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReportBuilder;

public class DiscreteSitesMatrixFELM extends AbstractFixedEffectLM {
	protected int numberOfSiteReportColumns;
	protected int numberOfAlleleReportColumns;
	protected double tol = 1e-10;
	protected ArrayList<DoubleMatrix> basisVectors;
	protected ArrayList<double[]> basisArrays;
	protected double[] normalizedData;
	protected int modeldf, errordf;
	protected boolean hasBasis = false;
	
	public DiscreteSitesMatrixFELM(Datum dataset, FixedEffectLMPlugin parentPlugin) {
		super(dataset, parentPlugin);
	}

	@Override
	public void initializeReportBuilders() {
		String tableName = "GLM Site Tests - " + myDatum.getName();
		String[] columnNames;
		if (permute) columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","perm_p","marker_df","error_df","model_df"};
		else columnNames = new String[]{"Trait","Marker","Chr","Position","marker_F","marker_p","marker_df","error_df","model_df"};
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
		if (!hasBasis) {
			hasBasis = true;
			calculateBasis();
		}
		byte[] geno = myGenoPheno.genotypeTable().genotypeAllTaxa(myCurrentSite);
		byte minor = myGenoPheno.genotypeTable().minorAllele(myCurrentSite);
		double[] coding = ModelEffectUtils.getNumericCodingForAdditiveModel(geno, minor);
		double[] normalGeno = normalizeVector(coding, true);
		
		for (int iter = 0; iter < 1; iter++) { //start for loop for testing

		double F,p;
		p = Double.NaN;
		int markerdf = 1;
		if (normalGeno.length == 0) {	//genotype df = 0, do not test
			F = 0;
			p = Double.NaN;
			markerdf = 0;
		} else {
			double corr = 0;
			for (int i = 0; i < normalGeno.length; i++) corr += normalGeno[i] * normalizedData[i];
			double corrsq = corr * corr;
			F = errordf * corrsq / (1 - corrsq );
//	        try {
//	        	p = LinearModelUtils.Ftest(F, 1, errordf);
//	        } catch (Exception e) {
//	        	p = Double.NaN;
//	        }
		}

		
		//report
        //add results to site report
        //columns = {"Trait","Marker","Chr","Position","marker_F","marker_p","perm_p","marker_df","error_df","model_df"}
		String siteName = myGenoPheno.genotypeTable().siteName(myCurrentSite);
        Object[] rowData = new Object[numberOfSiteReportColumns];
        int columnCount = 0;
        rowData[columnCount++] = currentTraitName;
        rowData[columnCount++] = siteName;
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomeName(myCurrentSite);
        rowData[columnCount++] = myGenoPheno.genotypeTable().chromosomalPosition(myCurrentSite);
        rowData[columnCount++] = new Double(F);
        rowData[columnCount++] = new Double(p);
        if (permute) rowData[columnCount++] = "";
        rowData[columnCount++] = new Integer(markerdf);
        rowData[columnCount++] = new Integer(errordf);
        rowData[columnCount++] = new Integer(modeldf - 1);
        siteReportBuilder.add(rowData);
        
		}//end for loop
	}

	public void calculateBasis() {
		//assume no missing data
		//create othornormal basis for the base model
		
		ArrayList<DoubleMatrix> modelVectors = new ArrayList<>();
		ArrayList<ModelEffect> effects = new ArrayList<>(myBaseModel);
		effects.remove(0); //do not use the mean
		for (ModelEffect effect : effects) {
			DoubleMatrix X = effect.getX();
			int ncol = X.numberOfColumns();
			int nrow = X.numberOfRows();
			for (int i = 0; i < ncol; i++) {
				DoubleMatrix aColumn = X.column(i);
				//center
				double mean = X.columnSum(i) / nrow;
				aColumn = aColumn.scalarAdd(-mean);
				modelVectors.add(aColumn);
			}
		}
		
		//for each vector v, subtract <u,v>*u for each of the previously calculated unit vectors
		basisVectors = new ArrayList<>();
		for (DoubleMatrix vector:modelVectors) {
			DoubleMatrix basis = vector.copy();
			for (DoubleMatrix u:basisVectors) {
				double alpha = u.crossproduct(vector).get(0, 0);
				basis.minusEquals(u.scalarMult(alpha));
			}
			double norm = Math.sqrt(basis.crossproduct(basis).get(0, 0));
			//if norm = 0, the vector is linearly dependent on previous vectors and is not part of the basis
			if (norm > tol) {
				basis.scalarMultEquals(1/norm);
				basisVectors.add(basis);
			}
		}

		basisArrays = new ArrayList<>();
		for (DoubleMatrix dm : basisVectors) {
			int n = dm.numberOfRows();
			double[] vals = new double[n];
			for (int i = 0; i < n; i++) vals[i] = dm.get(i, 0);
			basisArrays.add(vals);
		}
		
		normalizedData = normalizeVector(allData, true);
		modeldf = basisArrays.size() + 2; //plus one for the mean and one for the snp
		errordf = allData.length - modeldf;
	}
	
	private double[] normalizeVector(double[] input, boolean center) {
		int n = input.length;
		if (center) {
			double mean = 0;
			for (int i = 0; i < n; i++) mean += input[i];
			mean /= (double) n;
			double[] centered = new double[n];
			for (int i = 0; i < n; i++) {
				centered[i] = input[i] - mean;
			}
			input = centered;
		} 
		double[] output = Arrays.copyOf(input, n);
		
		for (double[] basis:basisArrays) {
			double alpha = 0;
			for (int i = 0; i < n; i++) alpha += input[i]*basis[i];
			for (int i = 0; i < n; i++) output[i] -= alpha*basis[i];
		}
		
		double norm = 0;
		for (int i = 0; i < n; i++) norm += output[i] * output[i];
		norm = Math.sqrt(norm);
		if (norm < tol) return new double[0];
		for (int i = 0; i < n; i++) output[i] /= norm;
		
		return output;
	}
	
	private double[] normalizeVector(float[] floatInput, boolean center) {
		
		int n = floatInput.length;
		double[] input = new double[n];
		if (center) {
			float mean = 0;
			for (int i = 0; i < n; i++) mean += floatInput[i];
			mean /= (double) n;
			for (int i = 0; i < n; i++) {
				input[i] = floatInput[i] - mean;
			}
		} else {
			for (int i = 0; i < n; i++) input[i] = floatInput[i];
		}
		
		double[] output = Arrays.copyOf(input, n);
		
		for (double[] basis:basisArrays) {
			double alpha = 0;
			for (int i = 0; i < n; i++) alpha += input[i]*basis[i];
			for (int i = 0; i < n; i++) output[i] -= alpha*basis[i];
		}
		
		double norm = 0;
		for (int i = 0; i < n; i++) norm += output[i] * output[i];
		norm = Math.sqrt(norm);
		if (norm < tol) return new double[0];
		for (int i = 0; i < n; i++) output[i] /= norm;
		
		return output;
	}

	@Override
	protected void getGenotypeAndUpdateMissing(BitSet missingObsBeforeSite) {
		missingObsForSite = new OpenBitSet(missingObsBeforeSite);
	}
	
}
