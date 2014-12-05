package net.maizegenetics.stats.linearmodels;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

public class SolveByOrthogonalizing {
	private List<ModelEffect> myBaseModel;
	private List<double[]> myBasisVectors;
	private List<double[]> myData;
	private List<double[]> myOrthogonalizedData;
	private SingularValueDecomposition baseSvd = null;
	private final static double tol = 1e-10;
	
	private SolveByOrthogonalizing() {
		
	}

	public static SolveByOrthogonalizing getInstanceFromModel(List<ModelEffect> baseModel, List<double[]> dataList) {
		SolveByOrthogonalizing sbo = new SolveByOrthogonalizing();
		sbo.myBaseModel = baseModel;
		sbo.myData = dataList;
		DoubleMatrix[][] design = sbo.createDesignMatricesFromModel();
		sbo.computeBaseSvd(design);
		sbo.OrthogonalizeData();
		return sbo;
	}
	
	public static SolveByOrthogonalizing getInstanceFromVectors(List<double[]> basisVectors, List<double[]> dataList) {
		SolveByOrthogonalizing sbo = new SolveByOrthogonalizing();
		sbo.myBasisVectors = basisVectors;
		sbo.myData = dataList;
		sbo.computeBaseSvd(sbo.createDesignMatricesFromVectors());
		sbo.OrthogonalizeData();
		return sbo;
	}
	
	public SolveByOrthogonalizing.Marker solveForR(SolveByOrthogonalizing.Marker marker) {
		if (marker.vector2 == null) return solveForR(marker.position(), marker.vector1());
		return solveForR(marker.position(), marker.vector1(), marker.vector2);
	}
	
	public SolveByOrthogonalizing.Marker solveForR(Position pos, double[] values) {
		double[] centeredValues = center(values);
		double[] orthoValues = orthogonalizeByBase(centeredValues);
		double[] csVales = centerAndScale(orthoValues);
		double[] orthogonalValues = centerAndScale(orthogonalizeByBase(center(values)));
		double[] rValues = myOrthogonalizedData.stream().mapToDouble(d -> innerProduct(d, orthogonalValues)).map(d -> d*d).toArray();
		
		int n = rValues.length;
		double[] pValues = new double[n]; 
		double modeldf = 1;
		if (baseSvd != null) modeldf += baseSvd.getRank();
		double errordf = values.length - modeldf - 1;

		for (int i = 0; i < n; i++) {
			pValues[i] = calculateP(calculateFfromR2(rValues[i], 1, errordf), 1, errordf);
		}
		return new SolveByOrthogonalizing.Marker(pos, rValues, pValues);
	}
	
	public SolveByOrthogonalizing.Marker solveForR(Position pos, double[] add, double[] dom) {
		if (dom == null) return solveForR(pos, add);
		double[] orthogonalAdd = orthogonalizeByBase(center(add));
		double[] orthogonalDom = orthogonalizeByBase(center(dom));
		
		//orthogonalize dom with respect to add
		double mult = innerProduct(orthogonalAdd, orthogonalDom)/innerProduct(orthogonalAdd, orthogonalAdd);
		int n = orthogonalDom.length;
		for (int i = 0; i < n; i++) {
			orthogonalDom[i] -= mult * orthogonalAdd[i];
		}
		
		//center and scale
		double[] v1 = centerAndScale(orthogonalAdd);
		double[] v2 = centerAndScale(orthogonalDom);
		double[] rValues = myOrthogonalizedData.stream()
				.mapToDouble(d -> {
					double r1 = innerProduct(v1, d);
					double r2 = innerProduct(v2, d);
					return r1 * r1 + r2 * r2;
				}).toArray();
		
		n = rValues.length;
		double[] pValues = new double[n]; 
		double modeldf = 2;
		if (baseSvd != null) modeldf += baseSvd.getRank();
		double errordf = v1.length - modeldf - 1;
		
		for (int i = 0; i < n; i++) {
			pValues[i] = calculateP(calculateFfromR2(rValues[i], 2, errordf), 2, errordf); 
		}
		return new SolveByOrthogonalizing.Marker(pos, rValues, pValues);
	}
	
	private DoubleMatrix[][] createDesignMatricesFromModel() {
		DoubleMatrix[][] designMatrices = new DoubleMatrix[1][];
		designMatrices[0] = myBaseModel.stream()
				.filter(a -> !a.getID().toString().toLowerCase().equals("mean"))
				.map(me -> me.getX())
				.toArray(DoubleMatrix[]::new);
		return designMatrices;
	}
	
	private DoubleMatrix[][] createDesignMatricesFromVectors() {
		DoubleMatrix[][] designMatrices = new DoubleMatrix[1][];
		designMatrices[0] = myBasisVectors.stream()
				.map(d -> DoubleMatrixFactory.DEFAULT.make(d.length, 1, d))
				.toArray(DoubleMatrix[]::new);
		return designMatrices;
	}
	
	private void computeBaseSvd(DoubleMatrix[][] designMatrices ) {
		if (designMatrices[0].length == 0 || designMatrices[0][0] == null) return; 
		
		DoubleMatrix X = DoubleMatrixFactory.DEFAULT.compose(designMatrices);
		
		//center the columns of X
		int nrows = X.numberOfRows();
		double dblnrows = nrows;
		int ncols = X.numberOfColumns();
		for (int c = 0; c < ncols; c++) {
			double mean = X.columnSum(c) / dblnrows;
			for (int r = 0; r < nrows; r++) {
				X.set(r, c, X.get(r, c) - mean);
			}
		}
		
		baseSvd = X.getSingularValueDecomposition();
	}
	
	private void OrthogonalizeData() {
		if (baseSvd == null) {
			myOrthogonalizedData = myData.stream()
					.map(d -> centerAndScale(Arrays.copyOf(d, d.length)))
					.collect(Collectors.toList());
			
		} else {
			myOrthogonalizedData = myData.stream()
					.map(d -> center(Arrays.copyOf(d, d.length)))
					.map(d -> orthogonalizeByBase(d))
					.map(d -> centerAndScale(d))
					.collect(Collectors.toList());
		}
	}
	
	private double[] orthogonalizeByBase(double[] vector) {
		if (baseSvd == null) {
			return Arrays.copyOf(vector, vector.length);
		}
		
		DoubleMatrix U = baseSvd.getU(false);
		int nrows = vector.length;
		double[] result = Arrays.copyOf(vector, nrows);
		int ncol = U.numberOfColumns();
		for (int i = 0; i < ncol; i++) {
			double[] u = U.column(i).to1DArray();
			double ip = innerProduct(vector, u);
			for (int j = 0; j < nrows; j++) result[j] -= ip * u[j];
		}
		return result;
	}
	
	public int baseDf() {
		int df = 1;  //the mean
		if (baseSvd != null) {
			double[] singularValues = baseSvd.getSingularValues();
			
			int n = singularValues.length;
			for (int i = 0; i < n; i++) {
				if (singularValues[i] > tol) df++;
			}
		}
		return df;
	}
	
	public static double innerProduct(double[] x, double[] y) {
		int n = x.length;
		return IntStream.range(0, n).mapToDouble(i -> x[i] * y[i]).sum();
	}
	
	public static double[] center(double[] values) {
		int n = values.length;
		double mean = Arrays.stream(values).sum() / n;
		for (int i = 0; i < n; i++) values[i] -= mean;
		return values;
	}
	
	public static double[] scale(double[] values) {
		int n = values.length;
		double divisor = Math.sqrt(innerProduct(values, values));
		for (int i = 0; i < n; i++) values[i] /= divisor;
		return values;
	}
	
	public static double[] centerAndScale(double[] values) {
		int n = values.length;
		double sum = 0;
		double sumsq = 0;
		for (int i = 0; i < n; i++) {
			double val = values[i];
			sum += val;
		}
		double mean = sum/n;
		
		for (int i = 0; i < n; i++) {
			values[i] = values[i] - mean;
			sumsq += values[i] * values[i];
		}
		
		double divisor = Math.sqrt(sumsq);
		for (int i = 0; i < n; i++) values[i] /= divisor;
		return values;
	}
	
	public static double calculateFfromR2(double r2, double markerDf, double errorDf) {
		return r2 / (1 - r2) * errorDf / markerDf;
	}
	
	public static double calculateP(double F, double markerDf, double errorDf) {
		if (!Double.isFinite(F)) return Double.NaN;
		double p;
 		try {
			p = LinearModelUtils.Ftest(F, markerDf, errorDf);
		} catch(Exception e) {
			p = Double.NaN;
		}
		return p;
	}
	
	public class Marker {
		Position myPosition;
		double[] vector1;
		double[] vector2;
		
		public Marker(Position pos, double[] values) {
			myPosition = pos;
			vector1 = values;
			vector2 = null;
		}
		
		public Marker(Position pos, double[] additive, double[] dominant) {
			myPosition = pos;
			vector1 = additive;
			vector2 = dominant;
		}
		
		public Position position() { return myPosition; }
		public double[] vector1() { return vector1; }
		public double[] vector2() { return vector2; }
	}
	
	
}
