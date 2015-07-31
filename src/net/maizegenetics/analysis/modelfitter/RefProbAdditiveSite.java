package net.maizegenetics.analysis.modelfitter;

import java.util.Arrays;

public class RefProbAdditiveSite extends AbstractAdditiveSite {

    private static final long serialVersionUID = 2040665024409852166L;
    private int ntaxa;
    private float[] cov;

    public RefProbAdditiveSite(int site, String chr, int pos, String id, CRITERION selectionCriteria, float[] covariate) {
        super(site, chr, pos, id, selectionCriteria);
        cov = covariate;
        ntaxa = cov.length;
    }

    @Override
    public double[] getCovariate() {
        double[] dcov = new double[ntaxa];
        for (int i = 0; i < ntaxa; i++)
            dcov[i] = cov[i];
        return dcov;
    }

    @Override
    public double[] getCovariate(int[] subset) {
        return Arrays.stream(subset).mapToDouble(i -> cov[i]).toArray();
    }

}
