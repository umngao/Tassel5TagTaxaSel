package net.maizegenetics.matrixalgebra.Matrix;

import java.util.Arrays;

import org.apache.log4j.Logger;

import net.maizegenetics.matrixalgebra.decomposition.BlasEigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.BlasSingularValueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.QRDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

public class BlasDoubleMatrix implements DoubleMatrix {
	private static Logger myLogger = Logger.getLogger(BlasDoubleMatrix.class);
	public static native void multMatrices(double[] A, int nrowsA, int ncolsA, double[] B, int ncolsB, double[] C, double alpha, double beta, 
			boolean transA, boolean transB); 
	public static native int solveLSdgelsd(double[] A, int Arows, int Acols, double[] B, int Bcols, double rcond, int[] rank);
	public static native int solveLSdgelsy(double[] A, int Arows, int Acols, double[] B, int Bcols, double rcond, int[] rank);
	public static native int singularValueDecompositionDgesvd(char jobu, char jobvt, int m, int n, double[] A, int lda, double[] S, double[] U, int ldu, double[] VT, int ldvt);
	public static native int singularValueDecompositionDgesdd(char jobz, int m, int n, double[] A, int lda, double[] S, double[] U, int ldu, double[] VT, int ldvt);
	public static native int eigenValueSymmetricDecomposition(int order, double[] A, double[] eigval); //A is the matrix on entry, eigenvectors on exit; returns error code
	
	//column major implementation
	protected double[] myMatrix;
	protected int nrows;
	protected int ncols;
	protected int size;
	
	static {
		System.load("/Users/pbradbury/Library/Developer/Xcode/DerivedData/MatrixAlgebraInterface-eeqitqwzxwogzifyqbgzjzcqincb/Build/Products/Debug/libMatrixAlgebraInterface.dylib");
	}
	
	public BlasDoubleMatrix() {
		
	}
	
	public BlasDoubleMatrix(int nrows, int ncols, double[] values, boolean columnMajor) {
		this.nrows = nrows;
		this.ncols = ncols;
		size = nrows * ncols;
		if (values == null) myMatrix = new double[size];
		else myMatrix = values;
		if (!columnMajor) transposeInPlace();
	}
	
	public BlasDoubleMatrix(double[][] values) {
		nrows = values[0].length;
		ncols = values.length;
		size = nrows * ncols;
		myMatrix = new double[size];
		int ptr = 0;
		for (int c = 0; c < ncols; c++) {
			for (int r = 0; r < nrows; r++) {
				myMatrix[ptr++] = values[c][r];
			}
		}
	}
	
	@Override
	public double get(int row, int col) {
		return myMatrix[getIndex(row,col)];
	}

	@Override
	public double getChecked(int row, int col) {
		if (row < nrows && col < ncols) {
			return myMatrix[getIndex(row,col)];
		}
		return Double.NaN;
	}

	@Override
	public void set(int row, int col, double value) {
		myMatrix[getIndex(row,col)] = value;
	}

	@Override
	public void setChecked(int row, int col, double value) {
		if (row < nrows && col < ncols) {
			myMatrix[getIndex(row,col)] = value;
		}
	}

	@Override
	public DoubleMatrix transpose() {
		BlasDoubleMatrix bdm = (BlasDoubleMatrix) copy();
		bdm.ncols = nrows;
		bdm.nrows = ncols;
		if (nrows > 1 && ncols > 1) {
			int ptr = 0;
			for (int i = 0; i < size; i++) {
				bdm.myMatrix[ptr] = myMatrix[i];
				ptr += bdm.nrows;
				if (ptr >= size) ptr -= size - 1;
			}
		}
		return bdm;
	}

	public void transposeInPlace() {
		BlasDoubleMatrix bdm = (BlasDoubleMatrix) transpose();
		nrows = bdm.ncols;
		ncols = bdm.nrows;
		myMatrix = bdm.myMatrix;
	}
	
	@Override
	public DoubleMatrix mult(DoubleMatrix dm, boolean transpose,
			boolean transposedm) {
		BlasDoubleMatrix B = (BlasDoubleMatrix) dm;
		BlasDoubleMatrix C = new BlasDoubleMatrix();
		
		if (transpose) C.nrows = ncols;
		else C.nrows = nrows;
		if (transposedm) C.ncols = B.nrows;
		else C.ncols = B.ncols;
		C.size = C.nrows * C.ncols;
		C.myMatrix = new double[C.size];
		
		multMatrices(myMatrix, nrows, ncols, B.myMatrix, B.nrows, C.myMatrix, 1.0, 0.0, transpose, transposedm);
		return C;
	}

	@Override
	public DoubleMatrix multadd(DoubleMatrix A, DoubleMatrix B, double alpha,
			double beta, boolean transpose, boolean transposeA) {
		BlasDoubleMatrix C = (BlasDoubleMatrix) A;
		BlasDoubleMatrix D = (BlasDoubleMatrix) B.copy();
		multMatrices(myMatrix, nrows, ncols, C.myMatrix, C.nrows, D.myMatrix, alpha, beta, transpose, transposeA);
		
		return D;
	}

	@Override
	public DoubleMatrix mult(DoubleMatrix dm) {
		return mult(dm, false, false);
	}

	@Override
	public DoubleMatrix crossproduct() {
		return mult(this, true, false);
	}

	@Override
	public DoubleMatrix crossproduct(DoubleMatrix dm) {
		return mult(dm, true, false);
	}

	@Override
	public DoubleMatrix tcrossproduct() {
		return mult(this, false, true);
	}

	@Override
	public DoubleMatrix tcrossproduct(DoubleMatrix dm) {
		return mult(dm, false, true);
	}

	@Override
	public DoubleMatrix concatenate(DoubleMatrix dm, boolean rows) {
		BlasDoubleMatrix B = (BlasDoubleMatrix) dm;
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		bdm.size = size + B.size;
		bdm.myMatrix = new double[bdm.size];
		
		if (rows & ncols == B.ncols) {
			bdm.nrows = nrows + B.nrows;
			bdm.ncols = ncols;
			int myptr = 0;
			int Bptr = 0;
			int bdmptr = 0;
			for (int c = 0; c < ncols; c++) {
				System.arraycopy(myMatrix, myptr, bdm.myMatrix, bdmptr, nrows);
				bdmptr += nrows;
				myptr += nrows;
				System.arraycopy(B.myMatrix, Bptr, bdm.myMatrix, bdmptr, B.nrows);
				bdmptr += B.nrows;
				Bptr += B.nrows;
			}
		
		} else if (nrows == B.nrows) {
			bdm.nrows = nrows;
			bdm.ncols = ncols + B.ncols;
			System.arraycopy(myMatrix, 0, bdm.myMatrix, 0, size);
			System.arraycopy(B.myMatrix, 0, bdm.myMatrix, size, B.size);
		}
		return null;
	}

	@Override
	public DoubleMatrix inverse() {
		BlasDoubleMatrix A = (BlasDoubleMatrix) this.copy();
		A.invert();
		return A;
	}

	@Override
	public void invert() {
		BlasDoubleMatrix inv = (BlasDoubleMatrix) generalizedInverseWithRank(new int[]{0});
		myMatrix = inv.myMatrix;
	}

	@Override
	public DoubleMatrix generalizedInverse() {
		return generalizedInverseWithRank(new int[]{0});
	}

	@Override
	public DoubleMatrix generalizedInverseWithRank(int[] rank) {
		//use solve
		return null;
	}

	@Override
	public DoubleMatrix solve(DoubleMatrix Y) {
		BlasDoubleMatrix bdy = (BlasDoubleMatrix) Y.copy();
		int[] rank = new int[]{0};
		int info = solveLSdgelsd(myMatrix, nrows, ncols, bdy.myMatrix, bdy.ncols, 1e-10, rank);
		if (info == 0) return bdy;
		myLogger.error(String.format("solve failed in BlasDoubleMatrix, info = %d\n",info));
		return null;
	}

	@Override
	public int numberOfRows() {
		return nrows;
	}

	@Override
	public int numberOfColumns() {
		return ncols;
	}

	@Override
	public DoubleMatrix row(int i) {
		if (i >= nrows) return null;
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		bdm.nrows = ncols;
		bdm.ncols = 1;
		bdm.size = ncols;
		bdm.myMatrix = new double[ncols];
		int myptr = i;
		int bdmptr = 0;
		while (myptr < size) {
			bdm.myMatrix[bdmptr++] = myMatrix[myptr];
			myptr += nrows;
		}
		return bdm;
	}

	@Override
	public DoubleMatrix column(int j) {
		if (j >= ncols) return null;
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		bdm.nrows = nrows;
		bdm.ncols = 1;
		bdm.size = nrows;
		int start = j*nrows;
		int end = start + nrows;
		bdm.myMatrix = Arrays.copyOfRange(myMatrix, start, end);
		return bdm;
	}

	@Override
	public DoubleMatrix[] getXtXGM() {
		DoubleMatrix xtx = crossproduct();
		DoubleMatrix g = inverse();
		BlasDoubleMatrix xg = (BlasDoubleMatrix) mult(g);
		BlasDoubleMatrix m = getIdentityMatrix(nrows);
		multMatrices(xg.myMatrix, xg.nrows, xg.ncols, myMatrix, nrows, m.myMatrix, -1, 1, false, true);
		
		return new DoubleMatrix[]{xtx, g, m};
	}

	@Override
	public DoubleMatrix copy() {
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		bdm.myMatrix = Arrays.copyOf(myMatrix, size);
		bdm.ncols = ncols;
		bdm.nrows = nrows;
		bdm.size = size;
		return bdm;
	}

	@Override
	public EigenvalueDecomposition getEigenvalueDecomposition() {
		return new BlasEigenvalueDecomposition(this);
	}

	@Override
	public SingularValueDecomposition getSingularValueDecomposition() {
		return new BlasSingularValueDecomposition(this);
	}

	@Override
	public QRDecomposition getQRDecomposition() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DoubleMatrix minus(DoubleMatrix dm) {
		BlasDoubleMatrix A = (BlasDoubleMatrix) copy();
		A.minusEquals(dm);
		return A;
	}

	@Override
	public void minusEquals(DoubleMatrix dm) {
		BlasDoubleMatrix bdm = (BlasDoubleMatrix) dm;
		for (int i = 0; i < size; i++) myMatrix[i] -= bdm.myMatrix[i];
	}

	@Override
	public DoubleMatrix plus(DoubleMatrix dm) {
		BlasDoubleMatrix A = (BlasDoubleMatrix) copy();
		A.plusEquals(dm);
		return A;
	}

	@Override
	public void plusEquals(DoubleMatrix dm) {
		BlasDoubleMatrix bdm = (BlasDoubleMatrix) dm;
		for (int i = 0; i < size; i++) myMatrix[i] += bdm.myMatrix[i];
	}

	@Override
	public DoubleMatrix scalarAdd(double s) {
		BlasDoubleMatrix A = (BlasDoubleMatrix) copy();
		A.scalarAddEquals(s);
		return A;
	}

	@Override
	public void scalarAddEquals(double s) {
		for (int i = 0; i < size; i++) myMatrix[i] += s;
	}

	@Override
	public DoubleMatrix scalarMult(double s) {
		BlasDoubleMatrix A = (BlasDoubleMatrix) copy();
		A.scalarMultEquals(s);
		return A;
	}

	@Override
	public void scalarMultEquals(double s) {
		for (int i = 0; i < size; i++) myMatrix[i] *= s;
	}

	@Override
	public DoubleMatrix getSelection(int[] rows, int[] columns) {
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		bdm.nrows = rows.length;
		bdm.ncols = columns.length;
		bdm.size = bdm.nrows * bdm.ncols;
		bdm.myMatrix = getSelectionFromDoubleArray(myMatrix, nrows, ncols, rows, columns);
		return bdm;
	}

	/**
	 * @param original	the original matrix as a one dimensional double array in column-major order
	 * @param nrows	the number of rows in the original matrix
	 * @param ncols	the number of columns in the original matrix
	 * @param rows	an int array of the selected rows
	 * @param columns	an int array of the selected columns
	 * @return the resulting double array, column-major order
	 */
	public static double[] getSelectionFromDoubleArray(double[] original, int nrows, int ncols, int[] rows, int[] columns) {
		if (rows == null && columns == null) return Arrays.copyOf(original, original.length);
		
		double[] selectedArray; 
		if (rows == null) {
			int numberOfSelectedRows = nrows;
			int numberOfSelectedColumns = columns.length;
			selectedArray = new double[numberOfSelectedRows * numberOfSelectedColumns];
			int ptr = 0;
			int myptr = 0;
			for (int c : columns) {
				int colstart = c * nrows;
				for (int r = 0; r < nrows; r++) {
					myptr = colstart + r;
					selectedArray[ptr++] = original[myptr];
				}
			}
		} else if (columns == null) {
			int numberOfSelectedRows = rows.length;
			int numberOfSelectedColumns = ncols;
			selectedArray = new double[numberOfSelectedRows * numberOfSelectedColumns];
			int ptr = 0;
			int myptr = 0;
			for (int c = 0; c < ncols; c++) {
				int colstart = c * nrows;
				for (int r : rows) {
					myptr = colstart + r;
					selectedArray[ptr++] = original[myptr];
				}
			}
		} else {
			int numberOfSelectedRows = rows.length;
			int numberOfSelectedColumns = columns.length;
			selectedArray = new double[numberOfSelectedRows * numberOfSelectedColumns];
			int ptr = 0;
			int myptr = 0;
			for (int c : columns) {
				int colstart = c * nrows;
				for (int r : rows) {
					myptr = colstart + r;
					selectedArray[ptr++] = original[myptr];
				}
			}
		}
		return selectedArray;
	}
	
	@Override
	public double rowSum(int row) {
		double sum = 0;
		int ptr = row;
		for (int c = 0; c < ncols; c++) {
			sum += myMatrix[ptr];
			ptr += nrows;
		}
		return sum;
	}

	@Override
	public double columnSum(int column) {
		double sum = 0;
		int start = column * nrows;
		int end = start + nrows;
		for (int ptr = start; ptr < end; ptr++) sum += myMatrix[ptr];
		return sum;
	}

	@Override
	public int columnRank() {
		BlasSingularValueDecomposition svd = new BlasSingularValueDecomposition(this, 'N');
		return svd.getRank();
	}

	public double[] getMatrix() { return myMatrix; }
	
	public double[] getMatrixCopy() { return Arrays.copyOf(myMatrix, size); }
	
	public int getSize() { return size; }
	
	public static BlasDoubleMatrix getDiagonalMatrix(double[] diag) {
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		int dim = diag.length;
		bdm.nrows = bdm.ncols = dim;
		bdm.size = bdm.nrows * bdm.ncols;
		bdm.myMatrix = new double[bdm.size];
		int ptr = 0;
		for (int i = 0; i < bdm.size; i += dim + 1)  bdm.myMatrix[i] = diag[ptr++];
		return bdm;
	}
	
	private int getIndex(int row, int col) {
		return col * nrows + row;
	}
	
	public static BlasDoubleMatrix getIdentityMatrix(int dim) {
		BlasDoubleMatrix bdm = new BlasDoubleMatrix();
		bdm.nrows = bdm.ncols = dim;
		bdm.size = dim * dim;
		bdm.myMatrix = new double[bdm.size];
		for (int i = 0; i < bdm.size; i += dim + 1)  bdm.myMatrix[i] = 1;
		return bdm;
	}
}
