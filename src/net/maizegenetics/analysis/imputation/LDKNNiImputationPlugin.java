package net.maizegenetics.analysis.imputation;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.MinMaxPriorityQueue;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.getDiploidValue;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.getUnphasedSortedDiploidValue;

/**
 * This imputation algorithm uses LD to identify good predictors for each SNP,
 * and then uses the high LD SNPs to identify K- Nearest Neighbors.
 * The genotype is called with a weighted mode of the KNNs.
 *
 * @author Daniel Money (Developer and developed the algorithm)
 * @author Ed Buckler (assisted in conversion to TASSEL)
 */
public class LDKNNiImputationPlugin extends AbstractPlugin {

    private PluginParameter<Integer> maxHighLDSites = new PluginParameter.Builder<>("maxHighLDSSites", 30, Integer.class)
            .range(Range.closed(2, 2000))
            .guiName("Max High LD Sites")
            .description("Maximum number of sites in high LD to use in imputation")
            .build();

    private PluginParameter<Integer> maxKNNTaxa = new PluginParameter.Builder<>("maxKNNTaxa", 10, Integer.class)
            .range(Range.closed(2, 200))
            .guiName("Max kNN taxa")
            .description("Maximum number of neighbours to use in imputation")
            .build();

    private PluginParameter<Integer> maxDistance = new PluginParameter.Builder<>("maxLDDistance", -1, Integer.class)
            .guiName("Max distance between site to find LD")
            .description("Maximum physical distance between sites to look for LD (-1 for no distance cutoff - unlinked chromosomes will be tested)")
            .build();


    private static final Logger myLogger = Logger.getLogger(LDKNNiImputationPlugin.class);

    public LDKNNiImputationPlugin() {
        super(null, false);
    }

    public LDKNNiImputationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() != 1) {
            throw new IllegalArgumentException("LDKNNiImputationPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        // Load in the genotype table
        Datum genoDatum=input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable)genoDatum.getData();

        // Create a multimap of the SNPs in highest LD with each SNP
/*Debatable on what to calc*/ //        Multimap<Position, Position> highLDMap = getHighLDMap(GenotypeTableBuilder.getHomozygousInstance(genotypeTable), maxHighLDSites());
        Multimap<Position, Position> highLDMap = getHighLDMap(genotypeTable, maxHighLDSites());
        System.out.println("LD calculated");

        GenotypeTableBuilder incSiteBuilder = GenotypeTableBuilder.getSiteIncremental(genotypeTable.taxa());
        //Start imputing site by site
        LongAdder sites1Kdone=new LongAdder();
        IntStream.range(0, genotypeTable.numberOfSites()).parallel().forEach(posIndex ->
        {
            Position position = genotypeTable.positions().get(posIndex);
            PositionList positionList = PositionListBuilder.getInstance(new ArrayList<>(highLDMap.get(position)));
            byte[] currGenos = genotypeTable.genotypeAllTaxa(posIndex);
            byte[] newGenos = new byte[currGenos.length];
            if (!genotypeTable.isPolymorphic(posIndex)) {
                byte monomorphicGenotype = getDiploidValue(genotypeTable.majorAllele(posIndex), genotypeTable.majorAllele(posIndex));
                for (int i = 0; i < newGenos.length; i++) {
                    newGenos[i] = (currGenos[i] == UNKNOWN_DIPLOID_ALLELE) ? monomorphicGenotype : currGenos[i];
                }
            } else {
                //create a site specific in memory filtered alignment to use
                final GenotypeTable ldGenoTable = GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(genotypeTable, positionList));
                for (int taxon = 0; taxon < currGenos.length; taxon++) {
                    if (currGenos[taxon] == UNKNOWN_DIPLOID_ALLELE) {
                        Multimap<Double, Byte> closeGenotypes = getClosestNonMissingTaxa(genotypeTable.taxa().get(taxon), genotypeTable,
                                ldGenoTable, position, maxKNNTaxa());
                        if (closeGenotypes.isEmpty()) {
                            newGenos[taxon] = currGenos[taxon];
                        } else {
                            newGenos[taxon] = impute(closeGenotypes, 20);
                        }
                    } else {
                        newGenos[taxon] = currGenos[taxon];
                    }

                }
            }
            incSiteBuilder.addSite(position, newGenos);
            if((posIndex+1)%100==0) {
                sites1Kdone.add(100);
                fireProgress(33+((int)(66*sites1Kdone.longValue())/genotypeTable.numberOfSites()));
            }
        });
        GenotypeTable impGenotypeTable = incSiteBuilder.build();

        return new DataSet(new Datum(genoDatum.getName()+"_KNNimp",impGenotypeTable,"Imputed genotypes by KNN imputation"),this);
    }

    private Multimap<Double, Byte> getClosestNonMissingTaxa(Taxon inputTaxon, GenotypeTable genotypeTable, GenotypeTable ldGenoTable,
                                                            Position targetPosition, int numberOfTaxa) {
        final int targetPosIdx = genotypeTable.positions().indexOf(targetPosition);
        final int inputTaxonIdx = genotypeTable.taxa().indexOf(inputTaxon);
        byte[] inputTaxonGenotypes = ldGenoTable.genotypeAllSites(inputTaxonIdx);
        MinMaxPriorityQueue<Tuple<Double, Byte>> topTaxa = IntStream.range(0, genotypeTable.numberOfTaxa())
                .filter(closeTaxonIdx -> closeTaxonIdx != inputTaxonIdx)
                .filter(closeTaxonIdx -> genotypeTable.genotype(closeTaxonIdx, targetPosIdx) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)
 /*Bad*/ //               .mapToObj(closeTaxonIdx -> new Tuple<>(IBSDistanceMatrix.computeHetDistances(inputTaxonGenotypes, ldGenoTable.genotypeAllSites(closeTaxonIdx), 10)[0], genotypeTable.genotype(closeTaxonIdx, targetPosIdx)))
 /*Best*/.mapToObj(closeTaxonIdx -> new Tuple<>(dist(inputTaxonGenotypes, ldGenoTable.genotypeAllSites(closeTaxonIdx), 10), genotypeTable.genotype(closeTaxonIdx, targetPosIdx)))
/*Best*/ //                .mapToObj(closeTaxonIdx -> new Tuple<>(LDKNNiImputationPluginEd.distance(inputTaxonGenotypes, ldGenoTable.genotypeAllSites(closeTaxonIdx), 10)[0], genotypeTable.genotype(closeTaxonIdx, targetPosIdx)))
                .filter(distanceTaxon -> !Double.isNaN(distanceTaxon.x))
                .collect(Collectors.toCollection(() -> MinMaxPriorityQueue.maximumSize(numberOfTaxa).create()));
        final Multimap<Double, Byte> distGenoMap = ArrayListMultimap.create();
        topTaxa.stream().forEach(distGeno -> distGenoMap.put(distGeno.x, distGeno.y));
        return distGenoMap;
    }

    private Multimap<Position, Position> getHighLDMap(GenotypeTable genotypeTable, int numberOfSNPs) {
        Multimap<Position, Position> highLDMap = ArrayListMultimap.create();
        final int numberOfSites = genotypeTable.numberOfSites();
        LongAdder sites1Kdone=new LongAdder();
        IntStream.range(0, genotypeTable.numberOfSites()).parallel()
                .forEach(posIndex ->
                {
                    MinMaxPriorityQueue<LDResult> highestLD = MinMaxPriorityQueue.orderedBy(LDResult.byR2Ordering.reverse())
                            .maximumSize(numberOfSNPs).create();
                    for (int site2 = 0; site2 < numberOfSites; site2++) {
                        if (posIndex == site2) {
                            continue;
                        }
                        if(maxDistance()>-1 && Math.abs(genotypeTable.chromosomalPosition(posIndex)-genotypeTable.chromosomalPosition(site2))>maxDistance()) {
                            continue;
                        }
                        LDResult ld = LinkageDisequilibrium.calculateBitLDForHaplotype(false, 20, genotypeTable, posIndex, site2);
                        if (Double.isNaN(ld.r2())) {
                            continue;
                        }
                        highestLD.add(ld);
                    }
                    List<Position> positionList = new ArrayList<>();
                    for (LDResult result : highestLD) {
                        positionList.add(genotypeTable.positions().get(result.site2()));
                    }
                    highLDMap.putAll(genotypeTable.positions().get(posIndex), positionList);
                    if((posIndex+1)%1000==0) {
                        sites1Kdone.add(1000);
                        fireProgress((int)(33*sites1Kdone.longValue())/numberOfSites);
                    }
                });

        return highLDMap;
    }

    /*
    Imputes to the most common genotype weighted by distance
     */
    static byte impute(Multimap<Double, Byte> distGeno, int useLDSites) {
        // Create an array to store the weighted counts
        double[] weightedCount = new double[256];

        // For each distance to genotype / genotype pair update the weighted counts
        distGeno.entries().forEach(entry -> {
            weightedCount[entry.getValue() + 128] += 1.0 / (1.0 + useLDSites * entry.getKey());
        });

        // Find the best genotype
        int bestGeno = 0;
        double bestWeightedCount = weightedCount[0];
        for (int i = 1; i < 256; i++) {
            if (weightedCount[i] > bestWeightedCount) {
                bestWeightedCount = weightedCount[i];
                bestGeno = i;
            }
        }
        return (byte) (bestGeno - 128);
    }

    /*
    Alternative approach to determine the best genotype
    This is an allele based approaches.  Provides a weighting of each allele, and Heterozygous is called currently
    if the second allele has the count more than 1/2 the most common.

     Provides roughly the same results as the genotype approach.  An ensemble is best but lower call rate.
     Both approaches have value.
     */
//    public static byte impute(Multimap<Double, Byte> distGeno, int useLDSites) {
//        // Create an array to store the weighted counts
//        TByteDoubleMap alleleWeight=new TByteDoubleHashMap(16);
//
//        // For each distance to genotype / genotype pair update the weighted counts
//        distGeno.entries().forEach(entry -> {
//            byte[] alleles= GenotypeTableUtils.getDiploidValues(entry.getValue());
//            double weight= 1.0 / (1.0 + useLDSites * entry.getKey());
//            alleleWeight.adjustOrPutValue(alleles[0],weight,weight);
//            alleleWeight.adjustOrPutValue(alleles[1],weight,weight);
//        });
//
//        double[] weights= alleleWeight.values();
//        Arrays.sort(weights);
//
//        if(weights.length==1) {
//            return GenotypeTableUtils.getUnphasedDiploidValue(alleleWeight.keys()[0], alleleWeight.keys()[0]);
//        } else if(weights[weights.length-1]>(2*weights[weights.length-2])) {
//            alleleWeight.retainEntries((a,w) -> w>=weights[weights.length-1]);
//            byte[] alleles=alleleWeight.keys();
//            return GenotypeTableUtils.getUnphasedDiploidValue(alleles[0], alleles[0]);
//        } else {
//            alleleWeight.retainEntries((a,w) -> w>=weights[weights.length-2]);
//            byte[] alleles=alleleWeight.keys();
//            return GenotypeTableUtils.getUnphasedDiploidValue(alleles[0], alleles[1]);
//        }
//    }

    @Override
    public String getCitation() {
        return "Daniel Money, Kyle Gardner, Heidi Schwaninger, Gan-Yuan Zhong, Sean Myles. (In Review) " +
                " LinkImpute: fast and accurate genotype imputation for non-model organisms";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "LD KNNi Imputation";
    }

    @Override
    public String getToolTipText() {
        return "LD KNNi Imputation";
    }


    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
//    public static void main(String[] args) {
//         GeneratePluginCode.generate(LDKNNiImputationPlugin.class);
//    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Maximum number of sites in high LD to use in imputation
     *
     * @return Max High LD Sites
     */
    public Integer maxHighLDSites() {
        return maxHighLDSites.value();
    }

    /**
     * Set Max High LD Sites. Maximum number of sites in high
     * LD to use in imputation
     *
     * @param value Max High LD Sites
     *
     * @return this plugin
     */
    public LDKNNiImputationPlugin maxHighLDSites(Integer value) {
        maxHighLDSites = new PluginParameter<>(maxHighLDSites, value);
        return this;
    }

    /**
     * Maximum number of neighbours to use in imputation
     *
     * @return Max kNN taxa
     */
    public Integer maxKNNTaxa() {
        return maxKNNTaxa.value();
    }

    /**
     * Set Max kNN taxa. Maximum number of neighbours to use
     * in imputation
     *
     * @param value Max kNN taxa
     *
     * @return this plugin
     */
    public LDKNNiImputationPlugin maxKNNTaxa(Integer value) {
        maxKNNTaxa = new PluginParameter<>(maxKNNTaxa, value);
        return this;
    }

    /**
     * Maximum physical distance between sites to look for
     * LD (-1 for no distance cutoff - unlinked chromosomes
     * will be tested)
     *
     * @return Max distance between site to find LD
     */
    public Integer maxDistance() {
        return maxDistance.value();
    }

    /**
     * Set Max distance between site to find LD. Maximum physical
     * distance between sites to look for LD (-1 for no distance
     * cutoff - unlinked chromosomes will be tested)
     *
     * @param value Max distance between site to find LD
     *
     * @return this plugin
     */
    public LDKNNiImputationPlugin maxDistance(Integer value) {
        maxDistance = new PluginParameter<>(maxDistance, value);
        return this;
    }


    /*
    Alternative to current IBS distance measure
    AA <> AA = 0
    Aa <> Aa = 0 distance (normal IBS distance this is 0.5)
    AA <> aa = 1 distance
     */
    //TODO Terry revisit where this should
    //TAS-787
    private static double dist(byte[] b1, byte[] b2, int min) {
        int d = 0;
        int c = 0;
        // Get the most similar snps to the current snp
        // Use the l most similar ones to calculate the distance
        for (int i = 0; i < b1.length; i++) {
            byte p1 = getUnphasedSortedDiploidValue(b1[i]);
            byte p2 = getUnphasedSortedDiploidValue(b2[i]);
            if ((p1 != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) && (p2 != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)) {
                // c counts how many snps we've actually used to scale the
                // distance with since some snps will be unknown
                c++;
                if (p1 != p2) {
                    if (GenotypeTableUtils.isHeterozygous(p1) ||
                            GenotypeTableUtils.isHeterozygous(p2)) {
                        d += 1;
                    } else {
                        d += 2;
                    }
                }
            }
        }
        // If across the l most similar snps there wasn't a single case
        // where both samples had a known genotype then set the distance to
        // NaN

        double ret;
        if (c < min) {
            ret = Double.NaN;
        }
        //Else return the scaled distance
        else {
            ret = ((double) d / (double) (2 * c));
        }

        return ret;
    }
}
