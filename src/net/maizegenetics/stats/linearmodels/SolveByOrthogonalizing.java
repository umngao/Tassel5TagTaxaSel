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
	List<double[]> myBasisVectors;
	List<double[]> myData;
	List<double[]> myOrthogonalizedData;
	SingularValueDecomposition baseSvd;
	
	private SolveByOrthogonalizing() {
		
	}

	public static SolveByOrthogonalizing getInstanceFromModel(List<ModelEffect> baseModel, List<double[]> dataList) {
		SolveByOrthogonalizing sbo = new SolveByOrthogonalizing();
		sbo.myBaseModel = baseModel;
		sbo.myData = dataList;
		sbo.computeBaseSvd(sbo.createDesignMatricesFromModel());
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
		return solveForR(marker.position(), marker.vector1());
	}
	
	public SolveByOrthogonalizing.Marker solveForR(Position pos, double[] values) {
		double[] centeredValues = center(values);
		double[] orthoValues = orthogonalizeByBase(centeredValues);
		double[] csVales = centerAndScale(orthoValues);
		double[] orthogonalValues = centerAndScale(orthogonalizeByBase(center(values)));
		double[] rValues = myOrthogonalizedData.stream().mapToDouble(d -> innerProduct(d, orthogonalValues)).toArray();
		return new SolveByOrthogonalizing.Marker(pos, rValues);
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
				.mapToDouble(d -> innerProduct(d, v1) + innerProduct(d, v2)).toArray();
		
		return new SolveByOrthogonalizing.Marker(pos, rValues);
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
				.map(d -> center(d))
				.map(d -> orthogonalizeByBase(d))
				.map(d -> centerAndScale(d))
				.collect(Collectors.toList());
	}
	
	private double[] orthogonalizeByBase(double[] vector) {
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
		double sum = 0;
		double sumsq = 0;
		for (int i = 0; i < n; i++) {
			double val = values[i];
			sum += val;
			sumsq += val * val;
		}
		double mean = sum/n;
		double divisor = Math.sqrt(sumsq);
		for (int i = 0; i < n; i++) values[i] = (values[i] - mean) / divisor;
		return values;
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
