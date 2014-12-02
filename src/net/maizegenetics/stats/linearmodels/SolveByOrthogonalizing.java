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
	List<ModelEffect> myBaseModel;
	List<double[]> myData;
	List<double[]> myOrthogonalizedData;
	SingularValueDecomposition baseSvd;
	
	public SolveByOrthogonalizing(List<ModelEffect> baseModel, List<double[]> dependentVariables) {
		myBaseModel = baseModel;
		myData = dependentVariables;
		computeBaseSvd();
		OrthogonalizeData();
	}
	
	public SolveByOrthogonalizing.Marker solveForR(SolveByOrthogonalizing.Marker marker) {
		return solveForR(marker.position(), marker.values());
	}
	
	public SolveByOrthogonalizing.Marker solveForR(Position pos, double[] values) {
		double[] orthogonalValues = centerAndScale(orthogonalizeByBase(values));
		double[] rValues = myOrthogonalizedData.stream().mapToDouble(d -> innerProduct(d, orthogonalValues)).toArray();
		return new SolveByOrthogonalizing.Marker(pos, rValues);
	}
	
	private void computeBaseSvd() {
		DoubleMatrix[][] designMatrices = new DoubleMatrix[1][];
		designMatrices[0] = (DoubleMatrix[]) myBaseModel.stream()
				.filter(a -> !a.getID().toString().toLowerCase().equals("mean"))
				.map(me -> me.getX())
				.toArray();
		
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
		myOrthogonalizedData = myData.stream()
				.map(d -> centerAndScale(orthogonalizeByBase(d)))
				.collect(Collectors.toList());
	}
	
	private double[] orthogonalizeByBase(double[] vector) {
		DoubleMatrix U = baseSvd.getU(false);
		int nrows = vector.length;
		double[] start = Arrays.copyOf(vector, nrows);
		return IntStream.range(0, U.numberOfColumns()).mapToObj(i -> U.column(i).to1DArray())
			.reduce(start, (a,b) -> {
				double ip = innerProduct(a, b);
				for (int j = 0; j < nrows; j++) a[j] = a[j] - b[j] * ip;
				return a;
			});
	}
	
	public static double innerProduct(double[] x, double[] y) {
		return IntStream.range(0, x.length).mapToDouble(i -> x[i] * y[i]).sum();
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
		double[] meanDiv = Arrays.stream(values).collect(() -> new double[]{0,0}, (a,b) -> {a[0] =+ b; a[1] += b * b;}, (a,b) -> {a[0] += b[0]; a[1] =+ b[1];});
		meanDiv[0] /= n;
		meanDiv[1] = Math.sqrt(meanDiv[1]);
		for (int i = 0; i < n; i++) values[i] = (values[i] - meanDiv[0]) / meanDiv[1];
		return values;
	}
	
	public class Marker {
		Position myPosition;
		double[] myValues;
		
		public Marker(Position pos, double[] values) {
			myPosition = pos;
			myValues = values;
		}
		
		public Position position() { return myPosition; }
		public double[] values() { return myValues; }
	}
	
	
}
