package net.maizegenetics.analysis.modelfitter;

import java.util.List;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;

public class GenotypeAdditiveSite extends AbstractAdditiveSite {

    private static final long serialVersionUID = -7891486608129027827L;
    private int[] bitStore;
    private double[] byteConversion;
    private int ntaxa;
    private int[] taxaIndex = null;

    /**
     * @param site						the site index from the originating GenotypeTable
     * @param genotype					genotypeAllTaxa for this site
     * @param majorAllele				major allele
     * @param majorAlleleFrequency		major allele frequency
     */
    public GenotypeAdditiveSite(int site, String chr, int pos, String id,
            CRITERION selectionCriterion, byte[] genotype,
            byte majorAllele, double majorAlleleFrequency) {
        super(site, chr, pos, id, selectionCriterion);
        byte unknown = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        ntaxa = genotype.length;
        double mean = 2 * majorAlleleFrequency;
        byteConversion = new double[] { -mean, 1 - mean, 2 - mean, 0 };

        int numberOfInts = ntaxa / 16;
        int remainder = ntaxa % 16;
        if (remainder > 0)
            numberOfInts++;

        bitStore = new int[numberOfInts];
        int intCount = 0;
        int genoCount = 0;

        while (genoCount < ntaxa) {
            int intStore = 0;
            for (int i = 0; i < 32; i += 2) {
                //calculate genoval
                int genoIndex;
                byte geno = genotype[genoCount++];
                if (geno == unknown) {
                    genoIndex = 3;
                } else {
                    genoIndex = 0;
                    byte[] alleles = GenotypeTableUtils.getDiploidValues(geno);
                    if (alleles[0] == majorAllele)
                        genoIndex++;
                    if (alleles[1] == majorAllele)
                        genoIndex++;
                }

                intStore = intStore | (genoIndex << i);
                if (genoCount == ntaxa)
                    break;
            }

            bitStore[intCount++] = intStore;
        }
    }

    @Override
    public void reindexTaxa(int[] taxaIndex, List<Integer> uniqueTaxa) {
        this.taxaIndex = taxaIndex;
        //count 0,1,2,3
        int[] counts = new int[4];
        for (Integer Ndx : uniqueTaxa) {
            counts[genotypeIndex(Ndx)]++;
        }

        double numerator = (double) (2 * counts[2] + counts[1]);
        double denominator = numerator + (counts[1] + 2 * counts[0]);
        double majorAlleleFreq = numerator / denominator;
        double mean = 2 * majorAlleleFreq;
        byteConversion = new double[] { -mean, 1 - mean, 2 - mean, 0 };
    }

    private int genotypeIndex(int n) {
        n = n * 2;
        int i = n / 32;
        int j = n % 32;
        int intStore = bitStore[i];
        return (intStore >> j) & 3;
    }

    @Override
    public double[] getCovariate() {
        if (taxaIndex == null)
            return getCovariateNoReindex();
        else
            return getCovariateWithReindex();
    }

    @Override
    public double[] getCovariate(int[] subset) {
        if (taxaIndex == null)
            return getCovariateNoReindex(subset);
        else
            return getCovariateWithReindex(subset);
    }

    public double[] getCovariateNoReindex() {
        double[] cov = new double[ntaxa];
        int intCount = 0;
        int genoCount = 0;
        while (genoCount < ntaxa) {
            int intStore = bitStore[intCount++];
            for (int i = 0; i < 32; i += 2) {
                int genoIndex = (intStore >> i) & 3;
                cov[genoCount++] = byteConversion[genoIndex];
                if (genoCount == ntaxa)
                    break;
            }
        }
        return cov;
    }

    public double[] getCovariateNoReindex(int[] subset) {
        int nobs = subset.length;
        double[] allCov = getCovariate();
        double[] cov = new double[nobs];
        for (int i = 0; i < nobs; i++)
            cov[i] = allCov[subset[i]];
        return cov;
    }

    public double[] getCovariateWithReindex() {
        int nval = taxaIndex.length;
        double[] cov = new double[nval];
        for (int i = 0; i < nval; i++) {
            cov[i] = byteConversion[genotypeIndex(taxaIndex[i])];
        }
        return cov;
    }

    public double[] getCovariateWithReindex(int[] subset) {
        int nval = subset.length;
        double[] cov = new double[nval];
        for (int i = 0; i < nval; i++) {
            cov[i] = byteConversion[genotypeIndex(taxaIndex[subset[i]])];
        }
        return cov;
    }
}
