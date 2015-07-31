package net.maizegenetics.analysis.modelfitter;

import java.io.Serializable;

public interface AdditiveSite extends Comparable<AdditiveSite>, Serializable {

    public static enum CRITERION {
        pval,
        aic,
        bic,
        mbic;
    }

    /**
     * @return		the covariate for this site
     */
    double[] getCovariate();

    /**
     * @param subset	an int array indexing a subset of taxa
     * @return			the covariate for the subset of taxa for this sites
     */
    double[] getCovariate(int[] subset);

    /**
     * @return	the site index corresponding to the site number in the source GenotypeTable
     */
    int siteNumber();

    /**
     * @return  the name of the chromosome of this site
     */
    String chromosomeName();
    
    /**
     * @return  the chromosomal position of this site
     */
    int position();
    
    /**
     * @return  the name (SNPID) of this site
     */
    String siteName();
    
    /**
     * @return     the value of the selection criterion for this site
     */
    double criterionValue();

    /**
     * @param value        the value of the selection criterion for this site
     */
    void criterionValue(double value);

    /**
     * @return     the selection criterion used for this site
     */
    CRITERION selectionCriterion();
}
