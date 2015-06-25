package net.maizegenetics.analysis.modelfitter;

import net.maizegenetics.analysis.modelfitter.AbstractForwardRegression.SiteInformation;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

public class AdditiveSite {
    private int siteIndex;
    private double SumSq;
    private int[] bitStore;
    private double[] byteConversion;
    private int ntaxa;
    
    public AdditiveSite() {
        siteIndex = -1;
        bitStore = null;
        byteConversion = null;
        SumSq = 0;
        ntaxa = 0;
    }
    
    public AdditiveSite(int site, byte[] genotype, byte majorAllele, double majorAlleleFrequency) {
        byte unknown = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        ntaxa = genotype.length;
        siteIndex = site;
        double mean = 2 * majorAlleleFrequency;
        byteConversion = new double[]{-mean, 1 - mean, 2 - mean, 0};
        
        int numberOfInts = ntaxa / 16;
        int remainder = ntaxa % 16;
        if (remainder > 0) numberOfInts++;
        
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
                    if (alleles[0] == majorAllele) genoIndex++;
                    if (alleles[1] == majorAllele) genoIndex++;
                }
                
                intStore = intStore | (genoIndex << i);
                if (genoCount == ntaxa) break;
            }
            
            bitStore[intCount++] = intStore;
        }
        
    }
    
    public double[] getCovariate() {
        double[] cov = new double[ntaxa];
        int intCount = 0;
        int genoCount = 0;
        while (genoCount < ntaxa) {
            int intStore = bitStore[intCount++];
            for (int i = 0; i < 32; i += 2) {
                int genoIndex = (intStore >> i) & 3;
                cov[genoCount++] = byteConversion[genoIndex];
                if (genoCount == ntaxa) break;
            }
        }
        return cov;
    }
    
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
                if (genoCount == nobs) break;
            }
        }
        return cov;
    }
    
    public int[] getIndex() {
        int[] ndx = new int[ntaxa];
        int intCount = 0;
        int genoCount = 0;
        while (genoCount < ntaxa) {
            int intStore = bitStore[intCount++];
            for (int i = 0; i < 32; i += 2) {
                int genoIndex = (intStore >> i) & 3;
                ndx[genoCount++] = genoIndex;
                if (genoCount == ntaxa) break;
            }
        }
        return ndx;
    }
    
    public int siteNumber() {
        return siteIndex;
    }
    
    public void SS(double ss) {
        SumSq = ss;
    }
    
    public double SS() {
        return SumSq;
    }
}
