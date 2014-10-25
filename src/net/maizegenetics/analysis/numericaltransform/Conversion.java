package net.maizegenetics.analysis.numericaltransform;

/**
 * User: dkroon
 * Date: Aug 1, 2005
 * Time: 10:41:05 AM
 */
public class Conversion {

    /**
     * Remove NaN values and resize the array.
     *
     * @param data Array containing NaN values.
     * @return Resized array without any NaN values.
     */
    public static double[] removeNaNValues(double[] data) {
        if (data == null || data.length == 0) {
            return null;
        }

        double[] cleanData = new double[data.length];
        int count = 0;
        for (int i = 0; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                cleanData[count++] = data[i];
            }
        }
        double[] tempData = new double[count];
        System.arraycopy(cleanData, 0, tempData, 0, count);
        return tempData;
    }

    /**
     * Remove NaN values and resize the array.
     *
     * @param data 2-D array containing NaN values.
     * @return Resized array without any NaN values.
     */
    public static double[][] removeNaNValues(double[][] data) {
        if (data == null || data.length == 0) {
            return null;
        }
        double[][] cleanData = new double[data.length][];
        int count = 0;
        for (int i = 0; i < data.length; i++) {
            cleanData[i] = removeNaNValues(data[i]);
        }

        return cleanData;
    }

    /**
     * Tests the normality of the given data and returns the pre-formatted results
     * from Kolmogorov-Smirnov, Cramer-vonMises and Anderson-Darling.
     * Utilizes the Stochastic Simulation in Java (SSJ) library.
     *
     * @param dataIn the data set to be tested
     * @return pre-formatted String of the normality test results
     */
    /*   public static String getNormalityResults(double[] data) {
    DoubleArrayList dal = new DoubleArrayList(data);

    ContinuousEmpiricalDist ced = new ContinuousEmpiricalDist(data);

    GofFormat.activeTests[2] = true;      // set Kolmogorov-Smirnov to be done
    GofFormat.activeTests[4] = true;      // set Cramer-vonMises to be done

    double[] results = new double[GofFormat.NTESTTYPES];
    double[] pResults = new double[GofFormat.NTESTTYPES];
    GofFormat.activeTests(dal, ced, results, pResults);

    return GofFormat.formatActiveTests(dal.size(), results, pResults);
    }
     */
    public static double[] normalizeData(double[] dataIn) {

        // translate into a 1-based array to prevent division by zero
        double[] data = new double[dataIn.length + 1];
        System.arraycopy(dataIn, 0, data, 1, dataIn.length);
        double avg = 0;
        double cumulativeValue = 0;
        int n = 0;
        for (int i = 1; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                cumulativeValue += data[i];
                n++;
            }
        }
        avg = cumulativeValue / n;
        double stDev = calculateStandardDeviation(dataIn);
        double[] result = new double[dataIn.length];
        for (int i = 0; i < dataIn.length; i++) {
            if (!Double.isNaN(dataIn[i])) {
                result[i] = (dataIn[i] - avg) / stDev;
            } else {
                result[i] = Double.NaN;
            }
        }
        return result;
    }

    /**
     * Normalizes the data in each column separately.
     *
     * @param dataIn
     * @return
     */
    public static double[][] normalizeData(double[][] dataIn) {
        double[][] result = new double[dataIn.length][dataIn[0].length];

        int colCount = dataIn[0].length;
        // do all of the calculations a column at a time
        for (int j = 0; j < colCount; j++) {
            double[] colData = new double[dataIn.length];
            for (int row = 0; row < dataIn.length; row++) {
                colData[row] = dataIn[row][j];
            }
            double[] normalizedData = normalizeData(colData);

            // fill the rows into the same column
            for (int q = 0; q < dataIn.length; q++) {
                if (!Double.isNaN(dataIn[q][j])) {
                    result[q][j] = normalizedData[q];
                } else {
                    result[q][j] = Double.NaN;
                }
            }
        }
        return result;
    }

    public static double calculateMean(double[] data) {
        if (data.length == 0) {
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                sum += data[i];
            }
        }
        return (double) sum / data.length;
    }

    /**
     * Calculate variance for a single array of values.
     *
     * @param dataIn
     * @return
     */
    public static double calculateVariance(double[] dataIn) {
        // translate into a 1-based array to prevent division by zero
        double[] data = new double[dataIn.length + 1];
        System.arraycopy(dataIn, 0, data, 1, dataIn.length);

        int n = 0;
        double sum = 0;
        double sumSqr = 0;

        for (int i = 1; i < data.length; i++) {

            if (!Double.isNaN(data[i])) {
                n += 1;
                sum += data[i];
                sumSqr += data[i] * data[i];
            }
        }

        return (sumSqr - sum * sum / n) / (n - 1);
    }

    /**
     *
     * @param data
     * @return
     */
    public static double calculateStandardDeviation(double[] data) {
        return Math.sqrt(calculateVariance(data));
    }
}
