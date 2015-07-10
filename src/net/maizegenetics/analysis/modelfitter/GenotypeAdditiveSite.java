package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;

public class GenotypeAdditiveSite extends AbstractAdditiveSite {
    private int[] bitStore;
    private double[] byteConversion;
    private int ntaxa;

    /**
     * @param site						the site index from the originating GenotypeTable
     * @param genotype					genotypeAllTaxa for this site
     * @param majorAllele				major allele
     * @param majorAlleleFrequency		major allele frequency
     */
    public GenotypeAdditiveSite(int site, CRITERION selectionCriterion, byte[] genotype,
            byte majorAllele, double majorAlleleFrequency) {
        super(site, selectionCriterion);
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
    public double[] getCovariate() {
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

    @Override
    public double[] getCovariate(int[] subset) {
        int nobs = subset.length;
        double[] cov = new double[nobs];
        int intCount = 0;
        int genoCount = 0;
        while (genoCount < nobs) {
            int intStore = bitStore[subset[intCount++]];
            for (int i = 0; i < 32; i += 2) {
                int genoIndex = (intStore >> i) & 3;
                cov[genoCount++] = byteConversion[genoIndex];
                if (genoCount == nobs)
                    break;
            }
        }
        return cov;
    }

}
