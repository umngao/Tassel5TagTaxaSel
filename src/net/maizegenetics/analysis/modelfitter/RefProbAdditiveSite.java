package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;

public class RefProbAdditiveSite extends AbstractAdditiveSite {
    private double SumSq;
    private double aic;
    private double bic;
    private double mbic;
    private double pval;
    private int ntaxa;
    private float[] cov;
    
    public RefProbAdditiveSite(int site, CRITERION selectionCriteria, float[] covariate) {
    	super(site, selectionCriteria);
    	cov = covariate;
    	ntaxa = cov.length;
    }
    
	@Override
	public double[] getCovariate() {
		double[] dcov = new double[ntaxa];
		for (int i = 0; i < ntaxa; i++) dcov[i] = cov[i];
		return dcov;
	}

	@Override
	public double[] getCovariate(int[] subset) {
		return Arrays.stream(subset).mapToDouble(i -> cov[i]).toArray();
	}



}
