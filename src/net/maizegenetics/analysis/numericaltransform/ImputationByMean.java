package net.maizegenetics.analysis.numericaltransform;

/**
 * Imputation of the missing data by mean of the respective column.
 *
 * @author Janu Verma
 *
 */
public class ImputationByMean {

    /**
     * Impute the missing values.
     *
     * @param data matrix
     * @return imputed data
     */
    public static double[][] impute(double[][] data) {
        int rows = data.length;
        int cols = data[0].length;
        //double[][] normData = Conversion.nomalizeData(data);
        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (!Double.isNaN(data[i][j])) {
                    result[i][j] = data[i][j];
                } else {
                    result[i][j] = mean(j, data);
                }
            }
        }
        return result;

    }

    /**
     *
     * @param col containing the missing data point.
     * @param data matrix
     * @return mean of the column.
     */
    public static double mean(int col, double[][] data) {
        int rows = data.length;
        double average;
        double sum = 0;
        for (int i = 0; i < rows; i++) {
            if (!Double.isNaN(data[i][col])) {
                sum = sum + data[i][col];
            }
        }
        average = sum / rows;
        return average;
    }

}
