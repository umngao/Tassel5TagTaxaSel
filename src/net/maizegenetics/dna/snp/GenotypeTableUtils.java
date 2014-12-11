/*
 * GenotypeTableUtils
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import java.nio.ByteBuffer;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import static net.maizegenetics.dna.snp.GenotypeTable.*;

/**
 * Utility methods for comparing, sorting, and counting genotypes.
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenotypeTableUtils {

    private static final Integer ONE = Integer.valueOf(1);
    private static final byte HIGHMASK = (byte) 0x0F;

    private GenotypeTableUtils() {
        // utility class
    }

    /**
     * This sorts alleles by frequency. Each cell in the given array contains a
     * diploid value which is separated and counted individually. Resulting
     * double dimension array holds alleles (bytes) in result[0]. And the counts
     * are in result[1]. Counts haploid values twice and diploid values once.
     * Higher ploids are not supported.
     *
     * @param data data
     *
     * @return alleles and counts
     */
    public static int[][] getAllelesSortedByFrequency(byte[] data) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < data.length; i++) {
            byte first = (byte) ((data[i] >>> 4) & 0xf);
            byte second = (byte) (data[i] & 0xf);
            if (first < RARE_ALLELE) {
                stateCnt[first]++;
            }
            if (second < RARE_ALLELE) {
                stateCnt[second]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    /**
     * This sorts alleles in a given site by frequency. Each cell in the given
     * array contains a diploid value which is separated and counted
     * individually. Resulting double dimension array holds alleles (bytes) in
     * result[0]. And the counts are in result[1]. Counts haploid values twice
     * and diploid values once. Higher ploids are not supported.
     *
     * @param data data
     * @param site site
     *
     * @return alleles and counts
     */
    public static int[][] getAllelesSortedByFrequency(byte[][] data, int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < data.length; i++) {
            byte first = (byte) ((data[i][site] >>> 4) & 0xf);
            byte second = (byte) (data[i][site] & 0xf);
            if (first < RARE_ALLELE) {
                stateCnt[first]++;
            }
            if (second < RARE_ALLELE) {
                stateCnt[second]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    public static byte[] getAlleles(byte[][] data, int site) {
        int[][] alleles = getAllelesSortedByFrequency(data, site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    public static Object[][] getAllelesSortedByFrequency(String[][] data, int site) {

        Map<String, Integer> stateCnt = new HashMap<>();
        for (int i = 0; i < data.length; i++) {
            String[] temp = data[i][site].split(":");
            String first;
            String second;
            if ((temp == null) || (temp.length == 0)) {
                first = second = UNKNOWN_ALLELE_STR;
            } else if (temp.length == 1) {
                first = second = temp[0].trim();
            } else {
                first = temp[0].trim();
                second = temp[1].trim();
            }
            if (!first.equalsIgnoreCase(UNKNOWN_ALLELE_STR)) {
                Integer count = stateCnt.get(first);
                if (count == null) {
                    stateCnt.put(first, ONE);
                } else {
                    stateCnt.put(first, count + 1);
                }
            }
            if (!second.equalsIgnoreCase(UNKNOWN_ALLELE_STR)) {
                Integer count = stateCnt.get(second);
                if (count == null) {
                    stateCnt.put(second, ONE);
                } else {
                    stateCnt.put(second, count + 1);
                }
            }
        }

        int count = stateCnt.size();

        Object[][] result = new Object[2][count];
        Iterator<String> itr = stateCnt.keySet().iterator();
        int index = 0;
        while (itr.hasNext()) {
            String key = itr.next();
            result[0][index] = key;
            result[1][index] = (Integer) stateCnt.get(key);
            index++;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if ((Integer) result[1][k] < (Integer) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    /**
     * Converts the byte representation of genotypes to string list of genotypes
     * @param data
     * @return list of the genotypes.
     */
    public static List<String> convertNucleotideGenotypesToStringList(byte[] data) {
        List<String> result=new ArrayList<>();
        for (byte b : data) {
            result.add(NucleotideAlignmentConstants.getHaplotypeNucleotide(b));
        }
        return result;
    }

    public static List<String> getAlleles(String[][] data, int site) {
        Object[][] alleles = getAllelesSortedByFrequency(data, site);
        if ((alleles == null) || (alleles.length == 0)) {
            return null;
        }
        int resultSize = alleles[0].length;
        List<String> result = new ArrayList<>();
        for (int i = 0; i < resultSize; i++) {
            result.add((String) alleles[0][i]);
        }
        return result;
    }

    public static int[][] getAllelesSortedByFrequency(GenotypeTable alignment, int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < alignment.numberOfTaxa(); i++) {
            byte[] dipB = alignment.genotypeArray(i, site);
            if (dipB[0] != GenotypeTable.UNKNOWN_ALLELE) {
                stateCnt[dipB[0]]++;
            }
            if (dipB[1] != GenotypeTable.UNKNOWN_ALLELE) {
                stateCnt[dipB[1]]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    public static Object[][] getDiploidsSortedByFrequency(GenotypeTable alignment, int site) {

        Integer ONE_INTEGER = 1;
        int numTaxa = alignment.numberOfTaxa();

        Map<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < numTaxa; r++) {
            String current = alignment.genotypeAsString(r, site);
            Integer num = diploidValueCounts.get(current);
            if (num == null) {
                diploidValueCounts.put(current, ONE_INTEGER);
            } else {
                diploidValueCounts.put(current, ++num);
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator<String> itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = itr.next();
            Integer count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

                if ((Integer) result[1][k] < (Integer) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    public static String[][] getAlleleStates(String[][] data, int maxNumAlleles) {

        int numSites = data[0].length;

        String[][] alleleStates = new String[numSites][16];
        for (int i = 0; i < numSites; i++) {
            for (int j = 0; j < 16; j++) {
                if (j == RARE_ALLELE) {
                    alleleStates[i][j] = GenotypeTable.RARE_ALLELE_STR;
                } else {
                    alleleStates[i][j] = UNKNOWN_ALLELE_STR;
                }
            }
        }

        for (int site = 0; site < numSites; site++) {
            List<String> alleles = GenotypeTableUtils.getAlleles(data, site);
            if (alleles != null) {
                int numAlleles = Math.min(alleles.size(), maxNumAlleles);
                for (int k = 0; k < numAlleles; k++) {
                    alleleStates[site][k] = (String) alleles.get(k);
                }
            }
        }

        return alleleStates;

    }

    public static byte[][] getDataBytes(String[] data) {

        int numTaxa = data.length;

        int numSites = data[0].length();

        byte[][] dataBytes = new byte[numTaxa][numSites];

        for (int site = 0; site < numSites; site++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                dataBytes[taxon][site] = NucleotideAlignmentConstants.getNucleotideDiploidByte(data[taxon].charAt(site));
            }
        }

        return dataBytes;

    }

    public static byte[][] getDataBytes(String[][] data, String[][] alleleStates, int maxNumAlleles) {

        int numTaxa = data.length;

        int numSites = data[0].length;

        byte[][] dataBytes = new byte[numTaxa][numSites];

        if (alleleStates.length == 1) {
            for (int site = 0; site < numSites; site++) {
                setDataBytes(data, alleleStates[0], maxNumAlleles, numTaxa, site, dataBytes);
            }
        } else {
            for (int site = 0; site < numSites; site++) {
                setDataBytes(data, alleleStates[site], maxNumAlleles, numTaxa, site, dataBytes);
            }
        }

        return dataBytes;

    }

    private static void setDataBytes(String[][] data, String[] alleleStates, int maxNumAlleles, int numTaxa, int site, byte[][] dataBytes) {
        if (data[0][0].contains(":")) {
            Pattern colon = Pattern.compile(":");
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                if (data[taxon][site].equalsIgnoreCase(UNKNOWN_DIPLOID_ALLELE_STR)) {
                    dataBytes[taxon][site] = UNKNOWN_DIPLOID_ALLELE;
                } else if (data[taxon][site].equals("?") || data[taxon][site].equals("?:?")) {
                    dataBytes[taxon][site] = UNKNOWN_DIPLOID_ALLELE;
                } else {
                    String[] siteval = colon.split(data[taxon][site]);
                    int[] byteval = new int[]{RARE_ALLELE, RARE_ALLELE};
                    for (int k = 0; k < maxNumAlleles; k++) {
                        if (alleleStates[k].equals(siteval[0])) {
                            byteval[0] = k;
                        }
                        if (alleleStates[k].equals(siteval[1])) {
                            byteval[1] = k;
                        }
                    }
                    dataBytes[taxon][site] = (byte) ((byteval[0] << 4) | byteval[1]);
                }
            }
        } else {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                if (data[taxon][site].equalsIgnoreCase(UNKNOWN_ALLELE_STR)) {
                    dataBytes[taxon][site] = UNKNOWN_DIPLOID_ALLELE;
                } else if (data[taxon][site].equals("?")) {
                    dataBytes[taxon][site] = UNKNOWN_DIPLOID_ALLELE;
                } else {
                    dataBytes[taxon][site] = RARE_DIPLOID_ALLELE;
                    for (byte k = 0; k < maxNumAlleles; k++) {
                        if (alleleStates[k].equals(data[taxon][site])) {
                            dataBytes[taxon][site] = (byte) (k | (k << 4));
                            break;
                        }
                    }
                }
            }
        }

    }

    /**
     * remove sites based on minimum frequency (the count of good bases,
     * INCLUDING GAPS) and based on the proportion of good alleles (including
     * gaps) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the
     * consensus
     * @param minimumCount minimum number of sequences with a good bases (not N
     * or ?), where GAP IS CONSIDERED A GOOD BASE
     */
    public static GenotypeTable removeSitesBasedOnFreqIgnoreMissing(GenotypeTable aa, double minimumProportion, double maximumProportion, int minimumCount) {
        int[] includeSites = getIncludedSitesBasedOnFreqIgnoreMissing(aa, minimumProportion, maximumProportion, minimumCount);
        GenotypeTable mlaa = FilterGenotypeTable.getInstance(aa, includeSites);
        return mlaa;
    }

    /**
     * get sites to be included based on minimum frequency (the count of good
     * bases, INCLUDING GAPS) and based on the proportion of good sites
     * (INCLUDING GAPS) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the
     * consensus
     * @param maximumProportion maximum proportion of sites different from the
     * consensus
     * @param minimumCount minimum number of sequences with a good base or a gap
     * (but not N or ?)
     */
    public static int[] getIncludedSitesBasedOnFreqIgnoreMissing(GenotypeTable aa, double minimumProportion, double maximumProportion, int minimumCount) {

        ArrayList<Integer> includeAL = new ArrayList<Integer>();
        for (int i = 0, n = aa.numberOfSites(); i < n; i++) {

            int totalNonMissing = aa.totalGametesNonMissingForSite(i);

            if ((totalNonMissing > 0) && (totalNonMissing >= (minimumCount * 2))) {

                double minorCount = aa.minorAlleleCount(i);
                double obsMinProp = 0.0;
                if (minorCount != 0) {
                    obsMinProp = minorCount / (double) totalNonMissing;
                }

                if ((obsMinProp >= minimumProportion) && (obsMinProp <= maximumProportion)) {
                    includeAL.add(i);
                }

            }
        }

        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        return includeSites;
    }

    /**
     * Remove sites based on site position (excluded sites are <firstSite and
     * >lastSite) This not effect any prior exclusions.
     *
     * @param aa the AnnotatedAlignment to filter
     * @param firstSite first site to keep in the range
     * @param lastSite last site to keep in the range
     */
    public static GenotypeTable removeSitesOutsideRange(GenotypeTable aa, int firstSite, int lastSite) {
        if ((firstSite < 0) || (firstSite > lastSite)) {
            return null;
        }
        if (lastSite > aa.numberOfSites() - 1) {
            return null;
        }
        return FilterGenotypeTable.getInstance(aa, firstSite, lastSite);
    }

    /**
     * Returns whether diploid allele values are heterozygous. First 4 bits in
     * byte is one allele value. Second 4 bits is other allele value.
     *
     * @param diploidAllele alleles
     *
     * @return true if allele values different; false if values the same.
     */
    public static boolean isHeterozygous(byte diploidAllele) {
        if (((diploidAllele >>> 4) & 0xf) == (diploidAllele & 0xf)) {
            return false;
        } else {
            return true;
        }
    }

    /**
     * Returns whether two diploid allele values are equal ignoring order.
     *
     * @param alleles1 diploid alleles 1
     * @param alleles2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqual(byte[] alleles1, byte[] alleles2) {

        if (((alleles1[0] == alleles2[0]) && (alleles1[1] == alleles2[1]))
                || ((alleles1[0] == alleles2[1]) && (alleles1[1] == alleles2[0]))) {
            return true;
        } else {
            return false;
        }

    }

    /**
     * Returns whether two diploid allele values are equal ignoring order.
     *
     * @param diploidAllele1 diploid alleles 1
     * @param diploidAllele2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqual(byte diploidAllele1, byte diploidAllele2) {

        if (diploidAllele1 != diploidAllele2) {
            byte reversed = (byte) ((diploidAllele1 << 4) | (diploidAllele1 >>> 4));
            if (reversed != diploidAllele2) {
                return false;
            }
        }

        return true;

    }

    /**
     * Returns whether two diploid allele values are equal ignoring order where
     * unknown values equal anything.
     *
     * @param alleles1 diploid alleles 1
     * @param alleles2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqualOrUnknown(byte[] alleles1, byte[] alleles2) {

        if (((alleles1[0] == GenotypeTable.UNKNOWN_ALLELE) && (alleles1[1] == GenotypeTable.UNKNOWN_ALLELE))
                || ((alleles2[0] == GenotypeTable.UNKNOWN_ALLELE) && (alleles2[1] == GenotypeTable.UNKNOWN_ALLELE))) {
            return true;
        }

        if (((alleles1[0] == alleles2[0]) && (alleles1[1] == alleles2[1]))
                || ((alleles1[0] == alleles2[1]) && (alleles1[1] == alleles2[0]))) {
            return true;
        } else {
            return false;
        }

    }

    /**
     * Returns whether two diploid allele values are equal ignoring order where
     * unknown values equal anything.
     *
     * @param diploidAllele1 diploid alleles 1
     * @param diploidAllele2 diploid alleles 2
     *
     * @return true if equal
     */
    public static boolean isEqualOrUnknown(byte diploidAllele1, byte diploidAllele2) {

        if ((diploidAllele1 == UNKNOWN_DIPLOID_ALLELE) || (diploidAllele2 == UNKNOWN_DIPLOID_ALLELE)) {
            return true;
        }

        if (diploidAllele1 != diploidAllele2) {
            byte reversed = (byte) ((diploidAllele1 << 4) | (diploidAllele1 >>> 4));
            if (reversed != diploidAllele2) {
                return false;
            }
        }

        return true;

    }

    /**
     * Return true if either at least one allele agree
     *
     * @param genotype1
     * @param genotype2
     * @return true if at least one allele is equal
     */
    public static boolean isPartiallyEqual(byte genotype1, byte genotype2) {
        int low1 = 0xF & genotype1;
        int low2 = 0xF & genotype2;
        if (low1 == low2) {
            return true;
        }
        int high1 = genotype1 >>> 4;
        if (high1 == low2) {
            return true;
        }
        int high2 = genotype2 >>> 4;
        if (low1 == high2) {
            return true;
        }
        if (high1 == high2) {
            return true;
        }
        return false;
    }

    public static boolean areEncodingsEqual(String[][][] encodings) {
        int numEncodings = encodings.length;
        for (int i = 1; i < numEncodings; i++) {
            int numSites = encodings[0].length;
            if (numSites != encodings[i].length) {
                return false;
            }
            for (int s = 0; s < numSites; s++) {
                int numCodes = encodings[0][s].length;
                if (numCodes != encodings[i][s].length) {
                    return false;
                }
                for (int c = 0; c < numCodes; c++) {
                    if (!encodings[0][s][c].equals(encodings[i][s][c])) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     * Combines two allele values into one diploid value. Assumed phased.
     *
     * @param a allele 1
     * @param b allele 2
     *
     * @return diploid value
     */
    public static byte getDiploidValuePhased(byte a, byte b) {
        return (byte) ((a << 4) | (HIGHMASK & b));
    }

    /**
     * Combines two allele values into one diploid value. Assumed phased.
     *
     * @param a allele 1
     * @param b allele 2
     *
     * @return diploid value
     */
    public static byte getDiploidValue(byte a, byte b) {
        return getDiploidValuePhased(a, b);
    }

    /**
     * Combines two allele values into one diploid value. In alphabetical order
     *
     * @param a allele 1
     * @param b allele 2
     *
     * @return diploid value sorted by order A < C < G < T
     */
    public static byte getUnphasedDiploidValue(byte a, byte b) {
        a = (byte) (HIGHMASK & a);
        b = (byte) (HIGHMASK & b);
        if (a < b) {
            return (byte) ((a << 4) | b);
        }
        return (byte) ((b << 4) | a);
    }

    /**
     * Combines two genotype values into one diploid value. Returns unknown if
     * either parent is heterozygous or unknown, or alleles are swapped.
     *
     * @param g1 genotype 1
     * @param g2 genotype 2
     * @return diploid value
     */
    public static byte getUnphasedDiploidValueNoHets(byte g1, byte g2) {
        if ((g2 == g1)&&(!isHeterozygous(g1))) {
            return g1;
        }
        if (g1 == UNKNOWN_DIPLOID_ALLELE) {
            return UNKNOWN_DIPLOID_ALLELE;
        }
        if (g2 == UNKNOWN_DIPLOID_ALLELE) {
            return UNKNOWN_DIPLOID_ALLELE;
        }
        if (isHeterozygous(g1)) {
            return UNKNOWN_DIPLOID_ALLELE;
        }
        if (isHeterozygous(g2)) {
            return UNKNOWN_DIPLOID_ALLELE;
        }

        return getUnphasedDiploidValue(g1, g2);
    }

    /**
     * Separates diploid allele value into it's two values.
     *
     * @param genotype diploid value
     *
     * @return separated allele values
     */
    public static byte[] getDiploidValues(byte genotype) {
        byte[] result = new byte[2];
        result[0] = (byte) ((genotype >>> 4) & 0xf);
        result[1] = (byte) (genotype & 0xf);
        return result;
    }

    /**
     * Method for getting TBits rapidly from major and minor allele arrays
     *
     * @param genotype
     * @param mjA
     * @param mnA
     * @return
     */
    public static BitSet[] calcBitPresenceFromGenotype(byte[] genotype, byte[] mjA, byte[] mnA) {
        int sites = genotype.length;
        if ((genotype.length != mjA.length) || (genotype.length != mnA.length)) {
            throw new ArrayIndexOutOfBoundsException("Input genotypes unequal in length");
        }
        OpenBitSet rMj = new OpenBitSet(genotype.length);
        OpenBitSet rMn = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            byte g = genotype[i];
            byte mj = mjA[i];
            byte mn = mnA[i];
            //           System.out.printf("inc:%d g:%d mj:%d mn:%d %n", i, g, mj, mn);
            if (mj == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == getDiploidValuePhased(mj, mj)) {
                rMj.fastSet(i);
                continue;
            }
            if (mn == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == getDiploidValuePhased(mn, mn)) {
                rMn.fastSet(i);
                continue;
            }
            byte het = getUnphasedDiploidValue(mj, mn);
            if (isEqual(g, het)) {
                rMj.fastSet(i);
                rMn.fastSet(i);
            }
        }
        return new BitSet[]{rMj, rMn};
    }
    
    public static BitSet calcBitPresenceFromGenotype(byte[] genotype, byte[] referenceValues) {
        int sites = genotype.length;
        if (genotype.length != referenceValues.length) {
            throw new ArrayIndexOutOfBoundsException("GenotypeTableUtils: calcBitPresenceFromGenotype: Input genotypes unequal in length");
        }
        OpenBitSet result = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            
            if (referenceValues[i] == UNKNOWN_ALLELE) {
                // do nothing
            } else if (referenceValues[i] == (byte) (genotype[i] & 0xf)) {
                result.fastSet(i);
            } else if (referenceValues[i] == (byte) ((genotype[i] >>> 4) & 0xf)) {
                result.fastSet(i);
            }

        }
        return result;
    }

    /**
     * Method for getting TBits rapidly from major and minor allele arrays
     *
     * @param genotype
     * @param mjA
     * @param mnA
     * @return
     */
    public static BitSet[] calcBitPresenceFromGenotype15(byte[] genotype, byte[] mjA, byte[] mnA) {
        int sites = genotype.length;
        if ((genotype.length != mjA.length) || (genotype.length != mnA.length)) {
            throw new ArrayIndexOutOfBoundsException("Input genotypes unequal in length");
        }
        ByteBuffer gBB = ByteBuffer.wrap(genotype);  //byte buffer make the code 20% faster for short sequences
        ByteBuffer mjBB = ByteBuffer.wrap(mjA);
        ByteBuffer mnBB = ByteBuffer.wrap(mnA);
        OpenBitSet rMj = new OpenBitSet(genotype.length);
        OpenBitSet rMn = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            byte g = gBB.get();
            byte mj = mjBB.get();
            byte mn = mnBB.get();
            // System.out.printf("inc:%d g:%d mj:%d mn:%d %n", i, g, mj, mn);
            if (mj == GenotypeTable.UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mj, mj)) {
                rMj.fastSet(i);
                continue;
            }
            if (mn == GenotypeTable.UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mn, mn)) {
                rMn.fastSet(i);
                continue;
            }
            byte het = GenotypeTableUtils.getUnphasedDiploidValue(mj, mn);
            if (GenotypeTableUtils.isEqual(g, het)) {
                rMj.fastSet(i);
                rMn.fastSet(i);
            }
        }
        return new BitSet[]{rMj, rMn};
    }

    /**
     * Method for getting Site Bits rapidly from major and minor alleles
     *
     * @param genotype
     * @param mj
     * @param mn
     * @return
     */
    public static BitSet[] calcBitPresenceFromGenotype(byte[] genotype, byte mj, byte mn) {
        int sites = genotype.length;
        OpenBitSet rMj = new OpenBitSet(genotype.length);
        OpenBitSet rMn = new OpenBitSet(genotype.length);
        for (int i = 0; i < sites; i++) {
            byte g = genotype[i];
            if (mj == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mj, mj)) {
                rMj.fastSet(i);
                continue;
            }
            if (mn == UNKNOWN_ALLELE) {
                continue;
            }
            if (g == GenotypeTableUtils.getDiploidValuePhased(mn, mn)) {
                rMn.fastSet(i);
                continue;
            }
            byte het = GenotypeTableUtils.getUnphasedDiploidValue(mj, mn);
            if (GenotypeTableUtils.isEqual(g, het)) {
                rMj.fastSet(i);
                rMn.fastSet(i);
            }
        }
        return new BitSet[]{rMj, rMn};
    }
    
    public static BitSet calcBitPresenceFromGenotype(byte[] genotype, byte referenceValue) {
        int sites = genotype.length;
        OpenBitSet result = new OpenBitSet(genotype.length);
        if (referenceValue == UNKNOWN_ALLELE) {
            return result;
        }
        for (int i = 0; i < sites; i++) {

            if (referenceValue == (byte) (genotype[i] & 0xf)) {
                result.fastSet(i);
            } else if (referenceValue == (byte) ((genotype[i] >>> 4) & 0xf)) {
                result.fastSet(i);
            }

        }
        return result;
    }
    
    public static float[][] convertGenotypeToFloatProbability(GenotypeTable genotype, boolean sitesByTaxa) {
    	int nsites = genotype.numberOfSites();
    	int ntaxa = genotype.numberOfTaxa();
    	float[][] out = null;
    	if (sitesByTaxa) {
    		out = new float[nsites][ntaxa];
    		for (int s = 0; s < nsites; s++) {
    			byte major = genotype.majorAllele(s);
    			for (int t = 0; t < ntaxa; t++) {
    				byte geno = genotype.genotype(t, s);
    				if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) out[s][t] = Float.NaN;
    				byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
    				if (alleles[0] == major) out[s][t] += 0.5;
    				if (alleles[1] == major) out[s][t] += 0.5;
    			}
    		} 
    	}
    	if (!sitesByTaxa) {
    		out = new float[ntaxa][nsites];
    		for (int s = 0; s < nsites; s++) {
    			byte major = genotype.majorAllele(s);
    			for (int t = 0; t < ntaxa; t++) {
    				byte geno = genotype.genotype(t, s);
    				if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) out[t][s] = Float.NaN;
    				byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
    				if (alleles[0] == major) out[t][s] += 0.5;
    				if (alleles[1] == major) out[t][s] += 0.5;
    			}
    		} 
    	}
    	return out;
	
    }

    public static double[][] convertGenotypeToDoubleProbability(GenotypeTable genotype, boolean sitesByTaxa) {
    	int nsites = genotype.numberOfSites();
    	int ntaxa = genotype.numberOfTaxa();
    	double[][] out = null;
    	if (sitesByTaxa) {
    		out = new double[nsites][ntaxa];
    		for (int s = 0; s < nsites; s++) {
    			byte major = genotype.majorAllele(s);
    			for (int t = 0; t < ntaxa; t++) {
    				byte geno = genotype.genotype(t, s);
    				if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) out[s][t] = Double.NaN;
    				byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
    				if (alleles[0] == major) out[s][t] += 0.5;
    				if (alleles[1] == major) out[s][t] += 0.5;
    			}
    		} 
    	}
    	if (!sitesByTaxa) {
    		out = new double[ntaxa][nsites];
    		for (int s = 0; s < nsites; s++) {
    			byte major = genotype.majorAllele(s);
    			for (int t = 0; t < ntaxa; t++) {
    				byte geno = genotype.genotype(t, s);
    				if (geno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) out[t][s] = Double.NaN;
    				byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
    				if (alleles[0] == major) out[t][s] += 0.5;
    				if (alleles[1] == major) out[t][s] += 0.5;
    			}
    		} 
    	}
    	return out;
    }
    
}
