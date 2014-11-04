package net.maizegenetics.analysis.numericaltransform;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Imputation of the missing data by k-nearest neighbors.
 *
 * @author Janu Verma
 *
 */
public class kNearestNeighbors {

    /**
     * Impute all the missing values.
     *
     * @param data matrix
     * @param k - number of nearest neighbors to be computed.
     * @param isManhattan - if true Manhattan distance will be used.
     * @param isCosine - if true Cosine similarity will be used.
     * @return imputed data matrix
     */
    public static double[][] impute(double[][] data, int k, boolean isManhattan, boolean isCosine) {
        int rows = data.length;
        int cols = data[0].length;
        double[][] result = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (!Double.isNaN(data[i][j])) {
                    result[i][j] = data[i][j];
                } else {
                    result[i][j] = calcKNN(data, i, j, k, isManhattan, isCosine);
                }
            }
        }
        return result;
    }

    /**
     * Compute the fill-in value.
     *
     * @param data matrix
     * @param row - row containing the missing value
     * @param col - col containing the missing value
     * @param k - number of nearest neighbors.
     * @param isManhattan - if true Manhattan distance will be used.
     * @param isCosine - if true Cosine similarity will be used.
     * @return value to be filled in for the missing data point.
     */
    public static double calcKNN(double[][] data, int row, int col, int k, boolean isManhattan, boolean isCosine) {
        double[][] neighbors;

        if (isCosine) {
            neighbors = cosineSimRank(data, row, col, k);
        } else {
            neighbors = KNearestNeighbor(data, row, col, k, isManhattan);
        }

        double num = 0.0;

        for (int i = 0; i < neighbors.length; i++) {
            if (!Double.isNaN(neighbors[i][col])) {
                num += neighbors[i][col];
            }
        }
        return (num / (k - 1));
    }

    /**
     * Computes the cosine similarity of two vectors.
     *
     * @param data1 array
     * @param data2 array
     * @return angle between the two vectors.
     */
    public static double cosine(double data1[], double data2[]) {
        int l1 = data1.length;
        int l2 = data2.length;

        double result = 0.0;

        for (int i = 0; i < l1; i++) {
            if ((!Double.isNaN(data1[i])) && (!Double.isNaN(data2[i]))) {
                result += (data1[i] * data2[i]);
            }
        }
        double norm1 = 0.0;
        for (int i = 0; i < l1; i++) {
            if (!Double.isNaN(data1[i])) {
                norm1 += Math.pow(data1[i], 2);
            }
        }

        double norm2 = 0.0;
        for (int j = 0; j < l2; j++) {
            if (!Double.isNaN(data2[j])) {
                norm2 += Math.pow(data2[j], 2);
            }
        }

        double normProduct = norm1 * norm2;

        return result / normProduct;

    }

    /**
     * Rank all the rows according to their cosine similarity with the given
     * row.
     *
     * @param data matrix
     * @param row
     * @param col
     * @param k
     * @return matrix whose rows are ordered according to their similarity with
     * given row.
     */
    public static double[][] cosineSimRank(double data[][], int row, int col, int k) {
        int nRows = data.length;
        int nCols = data[0].length;

        double[] query = data[row];
        double[][] neighbors = new double[k][nCols];

        Map<Integer, Double> distances = new HashMap<Integer, Double>();

        for (int i = 0; i < nRows; i++) {
            if (i != col) {
                double[] alpha = data[i];
                distances.put(i, cosine(query, alpha));
            }
        }

        for (int n = 0; n < k; n++) {
            double longestDistance = -10000.0;
            int IndexS = 0;

            Iterator it = distances.entrySet().iterator();
            while (it.hasNext()) {
                Map.Entry<Integer, Double> distance = (Map.Entry<Integer, Double>) it.next();
                if (distance.getValue() > longestDistance) {
                    longestDistance = distance.getValue();
                    IndexS = distance.getKey();

                }
            }
            neighbors[n] = data[IndexS];
            distances.remove(IndexS);
        }
        return neighbors;
    }

    /**
     * Compute Manhattan or Eulcidean distance between two vectors.
     *
     * @param data1 array
     * @param data2 array
     * @param isManhattan - true for Manhattan distance choice.
     * @return value of the distance between two vectors.
     */
    public static double distance(double data1[], double data2[], boolean isManhattan) {
        int l1 = data1.length;
        int l2 = data2.length;

        double result = 0.0;
        double count = 0.0;

        for (int i = 0; i < l1; i++) {
            if ((!Double.isNaN(data1[i])) && (!Double.isNaN(data2[i]))) {
                count += 1.0;
                if (isManhattan) {
                    result += Math.abs(data1[i] - data2[i]);
                } else {
                    result += Math.pow((data1[i] - data2[i]), 2);
                }
            }
        }
        return result / count;
    }

    /**
     * Rank the rows based on their distance from the given row.
     *
     * @param data matrix
     * @param row
     * @param col
     * @param k
     * @param isManhattan
     * @return matrix
     */
    public static double[][] KNearestNeighbor(double data[][], int row, int col, int k, boolean isManhattan) {
        int nRows = data.length;
        int nCols = data[0].length;
        double[] query = data[row];
        double[][] neighbors = new double[k][nCols];

        Map<Integer, Double> distances = new HashMap<Integer, Double>();

        for (int i = 0; i < nRows; i++) {
            double[] alpha = data[i];
            distances.put(i, distance(query, alpha, isManhattan));
        }

        for (int n = 0; n < k; n++) {
            double shortestDistance = -1.0;
            int IndexS = 0;

            Iterator it = distances.entrySet().iterator();
            while (it.hasNext()) {
                Map.Entry<Integer, Double> distance = (Map.Entry<Integer, Double>) it.next();
                if (distance.getValue() < shortestDistance || shortestDistance < 0) {
                    shortestDistance = distance.getValue();
                    IndexS = distance.getKey();
                }
            }
            neighbors[n] = data[IndexS];
            distances.remove(IndexS);
        }
        return neighbors;
    }

}
