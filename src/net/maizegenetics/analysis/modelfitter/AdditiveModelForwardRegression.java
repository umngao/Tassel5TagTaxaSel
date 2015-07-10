package net.maizegenetics.analysis.modelfitter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.apache.commons.math3.distribution.FDistribution;

import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.PartitionedLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;

public class AdditiveModelForwardRegression extends AbstractForwardRegression {
    double highestSS;
    int bestSite;
    
    public AdditiveModelForwardRegression(GenotypePhenotype data, int phenotypeIndex, double enterLimit, int maxVariants) {
        super(data, phenotypeIndex, enterLimit, maxVariants);
    }

    @Override
    public void fitModel() {
        int maxModelSize = myModel.size() + maxVariants;
        while (forwardStepParallel(true) && myModel.size() < maxModelSize);
    }
    
    @Override
    public void fitModelForSubsample(int[] subSample) {
        int numberOfSamples = subSample.length;
        
        //create myModel from myBaseModel for this subsample
        myModel = myBaseModel.stream().map(me -> me.getSubSample(subSample)).collect(Collectors.toList());
        double[] original = y;
        y = Arrays.stream(subSample).mapToDouble(i -> original[subSample[i]]).toArray();
        
        int maxModelSize = myModel.size() + maxVariants;
        while (forwardStepParallel(subSample, true)  && myModel.size() < maxModelSize);
        
        y = original;
    }

    private boolean forwardStep() {
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        PartitionedLinearModel plm = new PartitionedLinearModel(myModel, sflm);
        highestSS = 0;
        bestSite = -1;
        
        for (int s = 0; s < numberOfSites; s++) testSiteAsCovariate(plm, s);
        plm.setModelSS(highestSS);
        double[] thisFp = plm.getFp();
        if (thisFp[1] <= enterLimit) {
            addVariant(bestSite, thisFp[1]);
            myModel.add(new CovariateModelEffect(covariateForSite(bestSite)));
            
            return true;
        } else {
            return false;
        }
    }

    private boolean forwardStepParallel(boolean doParallel) {
        
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        PartitionedLinearModel plm = new PartitionedLinearModel(myModel, sflm);
        
        AdditiveSite bestSite = StreamSupport.stream(new ForwardStepAdditiveSpliterator(siteList, myModel, y), doParallel)
                .max((a,b) -> a.compareTo(b)).get();
        
        ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate());
        myModel.add(siteEffect);
        sflm = new SweepFastLinearModel(myModel, y);
        double[] errorSSdf = sflm.getResidualSSdf();
        double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
        double F,p;
        if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
            F = Double.NaN;
            p = Double.NaN;
        } else {
            F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
            p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
        }
        
        if (!Double.isNaN(p) && p <= enterLimit) {
            myFittedVariants.add(new Object[]{traitname, myGenotype.positions().get(bestSite.siteNumber()), new Integer(bestSite.siteNumber()), new Double(p)});  //Position, index, p-value
            addVariant(bestSite.siteNumber(), p);
            return true;
        }
        else return false;
    }
    
    private boolean forwardStepParallel(int[] subset, boolean doParallel) {
        
        AdditiveSite bestSite = StreamSupport.stream(new ForwardStepSubsettingAdditiveSpliterator(siteList, myModel, y, subset), doParallel)
                .max((a,b) -> a.compareTo(b)).get();

        ModelEffect siteEffect = new CovariateModelEffect(bestSite.getCovariate());
        myModel.add(siteEffect);
        SweepFastLinearModel sflm = new SweepFastLinearModel(myModel, y);
        double[] errorSSdf = sflm.getResidualSSdf();
        double[] siteSSdf = sflm.getIncrementalSSdf(myModel.size() - 1);
        double F,p;
        if (siteSSdf[1] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY || errorSSdf[0] < FDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY) {
            F = Double.NaN;
            p = Double.NaN;
        } else {
            F = siteSSdf[0] / siteSSdf[1] / errorSSdf[0] * errorSSdf[1];
            p = 1 - (new FDistribution(siteSSdf[1], errorSSdf[1]).cumulativeProbability(F));
        }
        
        if (!Double.isNaN(p) && p <= enterLimit) {
            myFittedVariants.add(new Object[]{traitname, myGenotype.positions().get(bestSite.siteNumber()), new Integer(bestSite.siteNumber()), new Double(p)});  //Position, index, p-value
            addVariant(bestSite.siteNumber(), p);
            return true;
        }
        else return false;
        
    }
    
    private void testSiteAsCovariate(PartitionedLinearModel plm, int site) {
        double ss = plm.testNewModelEffect(covariateForSite(site));
        if (ss > highestSS) {
            highestSS = ss;
            bestSite = site;
        }
    }

}
