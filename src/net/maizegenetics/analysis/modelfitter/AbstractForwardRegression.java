package net.maizegenetics.analysis.modelfitter;

import jal.Object.BinaryOperator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.analysis.modelfitter.AdditiveSite.CRITERION;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.stats.linearmodels.PartitionedLinearModel;
import net.maizegenetics.stats.linearmodels.SweepFastLinearModel;
import net.maizegenetics.util.Tuple;

/**
 * @author pbradbury
 *
 */
public abstract class AbstractForwardRegression implements ForwardRegression {
    //performs forward regression
    //requires one phenotype as a double array and a genotype table
    //accepts additional fixed effects, either covariates or factors
    //returns a list of markers with p-values
    //no missing values allowed in phenotype, factors, or covariates

    protected double[] y;     //data for a single phenotype (no missing data allowed)
    protected final GenotypePhenotype myGenotypePhenotype;
    protected final GenotypeTable myGenotype;
    protected final Phenotype myPhenotype;
    protected final double enterLimit;
    protected final int maxVariants;
    protected final int numberOfSites;
    protected final int numberOfObservations;
    protected List<AdditiveSite> siteList;
    protected final List<ModelEffect> myBaseModel;
    protected List<ModelEffect> myModel;
    protected List<Object[]> myFittedVariants = new ArrayList<>();
    protected String traitname;

    /**
     * @param data              a GenotypePhenotype object
     * @param phenotypeIndex    the attribute index of the phenotype to be analyzed
     * @param baseModel         the fixed effects in the model. If null, will be set to the mean only.
     * @param enterLimit        terms will be added to the model as long as the p-value of the next term is less than or equal to enterLimit
     * @param maxVariants       at most maxVariant terms will be fit. If the enterLimit is reached first, fewer than maxVariant terms will be fit.
     */
    public AbstractForwardRegression(GenotypePhenotype data, int phenotypeIndex, double enterLimit, int maxVariants) {
        myGenotypePhenotype = data;
        myGenotype = myGenotypePhenotype.genotypeTable();
        myPhenotype = myGenotypePhenotype.phenotype();
        traitname = myPhenotype.attributeName(phenotypeIndex);
        numberOfSites = data.genotypeTable().numberOfSites();
        numberOfObservations = data.numberOfObservations();
        myBaseModel = getBaseModel();
        myModel = new ArrayList<>(myBaseModel);
        y = AssociationUtils.convertFloatArrayToDouble((float[]) data.phenotype().attribute(phenotypeIndex).allValues());
        
        this.enterLimit = enterLimit;
        this.maxVariants = maxVariants;
        
        //Initialize the siteList
        siteList = new ArrayList<>();
        for (int s = 0; s < numberOfSites; s++) 
            siteList.add(new GenotypeAdditiveSite(s, CRITERION.ss, myGenotype.genotypeAllTaxa(s), myGenotype.majorAllele(s), myGenotype.majorAlleleFrequency(s))); 
    }
    
    protected List<ModelEffect> getBaseModel() {
        List<ModelEffect> base = new ArrayList<>();
        int[] mean = new int[numberOfObservations];
        ModelEffect meanEffect = new FactorModelEffect(mean, false, "mean");
        base.add(meanEffect);
        for (PhenotypeAttribute factor : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor) ) {
            ModelEffect factorEffect = new FactorModelEffect( ((CategoricalAttribute) factor).allIntValues(), true, factor.name());
            base.add(factorEffect);
        }
        for (PhenotypeAttribute cov : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate) ) {
            double[] values = AssociationUtils.convertFloatArrayToDouble(((NumericAttribute) cov).floatValues());
            ModelEffect covEffect = new CovariateModelEffect(values, cov.name());
            base.add(covEffect);
        }
        return base;
    }
    
    @Override
    public List<Object[]> fittedModel() {
        return myFittedVariants;
    }

    protected void addVariant(int site, double p) {
        Position myPosition = myGenotype.positions().get(site);
        myFittedVariants.add(new Object[]{traitname, myPosition.getSNPID(), myPosition.getChromosome().getName(), myPosition.getPosition(), p, -Math.log10(p)});
    }
    
    /**
     * @param site      a site in the GenotypeTable
     * @return          a double array of counts of the major allele at this site (0,1,2). Missing values are set to the site mean (2 * major allele frequency).
     */
    public double[] covariateForSite(int site) {
        double[] siteCov = new double[numberOfObservations];
        byte[] geno = myGenotypePhenotype.genotypeAllTaxa(site);
        byte majorAllele =   myGenotype.majorAllele(site);
        double majorAlleleFreq = myGenotype.majorAlleleFrequency(site);
        byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        for (int t = 0; t < numberOfObservations; t++) {
            if (geno[t] == NN)
                siteCov[t] = 2 * majorAlleleFreq;
            else {
                byte[] alleles = GenotypeTableUtils.getDiploidValues(geno[t]);
                if (alleles[0] == majorAllele)
                    siteCov[t]++;
                if (alleles[1] == majorAllele)
                    siteCov[t]++;
            }
        }
        return siteCov;
    }

    public Stream<SiteInformation> streamCovariateForSite() {
        return IntStream.range(0, numberOfSites).mapToObj(s -> new SiteInformation(s, covariateForSite(s), 0));
    }
    
    class SiteInformation {
        int siteIndex;
        double[] covariate;
        double SumSq;
        
        SiteInformation() {
            siteIndex = -1;
            covariate = null;
            SumSq = 0;
        }
        
        SiteInformation(int site, double[] cov, double ss) {
            siteIndex = site;
            covariate = cov;
            SumSq = ss;
        }
        
        SiteInformation best(SiteInformation otherSite) {
            if (otherSite.SumSq > SumSq) return otherSite;
            return this;
        }
        
    }
    
    public class ForwardStepAdditiveSpliterator implements Spliterator<AdditiveSite> {
        private final PartitionedLinearModel plm;
        private List<AdditiveSite> mySites;
        private final List<ModelEffect> baseModel;
        private final double[] y;
        private int origin;
        private final int end;
        
        public ForwardStepAdditiveSpliterator(List<AdditiveSite> siteList, List<ModelEffect> baseModel, double[] y) {
            SweepFastLinearModel sflm = new SweepFastLinearModel(baseModel, y);
            plm = new PartitionedLinearModel(baseModel, sflm);
            mySites = siteList;
            this.baseModel = baseModel;
            this.y = y;
            origin = 0;
            end = siteList.size();
        }
        
        @Override
        public boolean tryAdvance(Consumer<? super AdditiveSite> action) {
            if (origin == end) return false;
            AdditiveSite as = mySites.get(origin);
            as.criterionValue(plm.testNewModelEffect(as.getCovariate()));
            action.accept(as);
            origin++;
            return true;
        }

        @Override
        public Spliterator<AdditiveSite> trySplit() {
            int numberRemaining = end - origin;
            if (numberRemaining < 5000) return null;
            int mid = origin + numberRemaining / 2;
            List<AdditiveSite> splitSublist = mySites.subList(origin, mid);
            origin = mid;
            return new ForwardStepAdditiveSpliterator(splitSublist, baseModel, y);
        }

        @Override
        public long estimateSize() {
            return end - origin;
        }

        @Override
        public int characteristics() {
            return Spliterator.IMMUTABLE + Spliterator.NONNULL + Spliterator.SIZED + Spliterator.SUBSIZED;
        }
        
    }
    
    public class ForwardStepAdditiveSubsettingSpliterator implements Spliterator<AdditiveSite> {
        private final PartitionedLinearModel plm;
        private List<AdditiveSite> mySites;
        private final List<ModelEffect> baseModel;
        private final double[] y;
        private int origin;
        private final int end;
        private final int[] subset;
        
        public ForwardStepAdditiveSubsettingSpliterator(List<AdditiveSite> siteList, List<ModelEffect> baseModel, double[] y, int[] subset) {
            SweepFastLinearModel sflm = new SweepFastLinearModel(baseModel, y);
            plm = new PartitionedLinearModel(baseModel, sflm);
            mySites = siteList;
            this.baseModel = baseModel;
            this.y = y;
            origin = 0;
            end = siteList.size();
            this.subset = subset;
        }
        
        @Override
        public boolean tryAdvance(Consumer<? super AdditiveSite> action) {
            if (origin == end) return false;
            AdditiveSite as = mySites.get(origin);
            as.criterionValue(plm.testNewModelEffect(as.getCovariate(subset)));
            action.accept(as);
            origin++;
            return true;
        }

        @Override
        public Spliterator<AdditiveSite> trySplit() {
            int numberRemaining = end - origin;
            if (numberRemaining < 5000) return null;
            int mid = origin + numberRemaining / 2;
            List<AdditiveSite> splitSublist = mySites.subList(origin, mid);
            origin = mid;
            return new ForwardStepAdditiveSubsettingSpliterator(splitSublist, baseModel, y, subset);
        }

        @Override
        public long estimateSize() {
            return end - origin;
        }

        @Override
        public int characteristics() {
            return Spliterator.IMMUTABLE + Spliterator.NONNULL + Spliterator.SIZED + Spliterator.SUBSIZED;
        }
        
    }
}
