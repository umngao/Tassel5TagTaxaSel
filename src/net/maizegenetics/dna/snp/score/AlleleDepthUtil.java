/*
 *  AlleleDepthUtil
 */
package net.maizegenetics.dna.snp.score;

import java.util.Arrays;

/**
 * Depth compression utility. Depths are scored in byte (256 states). 0 to 127
 * are used for exact depth values. -128 to -1 are used for log approximations
 * (1.0746^(-x)). This permits depths upto 10,482 to be stored with exact
 * precision for low values, and <3% error for high depths.
 *
 * @author Terry Casstevens
 * @author Robert Bukowski
 */
public class AlleleDepthUtil {

    private static final double LOG_BASE = 1.0746;  // LOG_BASE^128 = 10,482
    private static final double R_LOG_CONV = 1.0 / Math.log(LOG_BASE);
    private static final double LOG_CONV = 1.0 / R_LOG_CONV;
    private static final int[] BYTE_TO_INT = new int[256];
    private static final int MAX_ACC_DEPTH = 182;
    private static final int MIN_ACC_BYTE = 127 - MAX_ACC_DEPTH;
    private static final int OFFSET = 126;
    private static final double ADJ = 0.5;

    static {
        Arrays.fill(BYTE_TO_INT, -1);
        for (int i = 0; i < 256; i++) {
            BYTE_TO_INT[i] = decode((byte) i);
        }
    }

    private AlleleDepthUtil() {
        // utility
    }

    /**
     * Returns the depth combined depths in byte format. Converts both to
     * integers, adds, and reconverts back.
     *
     * @param depth1 first depth in byte format
     * @param depth2 second depth in byte format
     * @return byte version of combined depth
     */
    public static byte addByteDepths(byte depth1, byte depth2) {
        return depthIntToByte(depthByteToInt(depth1) + depthByteToInt(depth2));
    }

    /**
     * Converts integer depth to byte version. Positive return values are exact,
     * negative return values are log approximations
     */
    public static byte depthIntToByte(int depth) {
        byte bdepth;
        int itd;

        if (depth <= 127) {
            itd = depth;
        } else if (depth <= MAX_ACC_DEPTH) {
            itd = 127 - depth;
        } else {
            itd = (int) (-R_LOG_CONV * Math.log(depth - OFFSET));
            if (itd < -128) {
                itd = -128;
            }
        }
        bdepth = (byte) itd;

        return bdepth;
    }

    /**
     * Converts integer array of depth to byte array version. Positive return
     * values are exact, negative return values are log approximations
     */
    public static byte[] depthIntToByte(int[] depth) {
        byte[] result = new byte[depth.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = depthIntToByte(depth[i]);
        }
        return result;
    }

    /**
     * Converts a two dimensional integer array of depths to a 2D byte array
     * version. Positive return values are exact, negative return values are log
     * approximations
     */
    public static byte[][] depthIntToByte(int[][] depth) {
        byte[][] result = new byte[depth.length][depth[0].length];
        for (int i = 0; i < result.length; i++) {
            result[i] = depthIntToByte(depth[i]);
        }
        return result;
    }

    /**
     * Converts byte depth to integer version. Positive depths value are exact,
     * negative depths are log approximations
     */
    public static int depthByteToInt(byte depth) {
        return BYTE_TO_INT[depth & 0xFF];
    }

    /**
     * Converts byte arrays of depth to integer array version. Positive depths
     * value are exact, negative depths are log approximations
     */
    public static int[] depthByteToInt(byte[] depth) {
        int[] result = new int[depth.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = depthByteToInt(depth[i]);
        }
        return result;
    }

    private static int decode(byte bdepth) {
        int depth;
        if (bdepth >= 0) {
            depth = bdepth;
        } else if (bdepth >= MIN_ACC_BYTE) {
            depth = 127 - bdepth;
        } else {
            depth = OFFSET + (int) (Math.exp(-LOG_CONV * (bdepth - ADJ)));
        }

        return depth;
    }

}
