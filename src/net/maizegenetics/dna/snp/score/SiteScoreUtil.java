/*
 *  SiteScoreUtil
 */
package net.maizegenetics.dna.snp.score;

import java.util.Arrays;

/**
 * @author Terry Casstevens
 */
public class SiteScoreUtil {

    private static final float[] BYTE_TO_FLOAT = new float[256];

    static {
        Arrays.fill(BYTE_TO_FLOAT, -1);
        for (int i = 0; i < 256; i++) {
            BYTE_TO_FLOAT[i] = decode((byte) i);
        }
    }

    private SiteScoreUtil() {
        // utility
    }

    /**
     * Converts float value to byte version.
     */
    public static byte floatToBytePercentage(float value) {
        if ((value < 0.0) || (value > 1.0)) {
            throw new IllegalArgumentException("SiteScoreUtil: floatToBytePercentage: value must be between 0.0 and 1.0");
        }

        return (byte) Math.round(255.0f * value);
    }

    public static byte[] floatToBytePercentage(float[] values) {
        byte[] result = new byte[values.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = floatToBytePercentage(values[i]);
        }
        return result;
    }

    public static byte[][] floatToBytePercentage(float[][] values) {
        byte[][] result = new byte[values.length][values[0].length];
        for (int i = 0; i < result.length; i++) {
            result[i] = floatToBytePercentage(values[i]);
        }
        return result;
    }

    /**
     * Converts byte value to float percentage.
     */
    public static float byteToFloatPercentage(byte value) {
        return BYTE_TO_FLOAT[value & 0xFF];
    }

    public static float[] byteToFloatPercentage(byte[] values) {
        float[] result = new float[values.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = byteToFloatPercentage(values[i]);
        }
        return result;
    }

    private static float decode(byte value) {
        return value / 255;
    }

}
