package net.maizegenetics.analysis.imputation;

import com.google.common.collect.*;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by edbuckler on 5/4/15.
 */
public class LDKNNiImputationPluginEd extends AbstractPlugin{
    private PluginParameter<String> hmpFile= new PluginParameter.Builder<>("i",null,String.class).guiName("Target file").inFile().required(true)
            .description("Input HapMap file of target genotypes to impute. Accepts all file types supported by TASSEL5.").build();
    private PluginParameter<String> outFileBase= new PluginParameter.Builder<>("o",null,String.class).guiName("Output filename").outFile().required(true)
            .description("Output file; hmp.txt.gz and .hmp.h5 accepted.").build();
    private PluginParameter<Integer> maxHighLDSites= new PluginParameter.Builder<>("maxHighLDSites",20,Integer.class).range(Range.closed(2,2000)).guiName("Maximum number of high LD sites")
            .description("Preferred haplotype block size in sites (use same as in FILLINFindHaplotypesPlugin)").build();
    private PluginParameter<Integer> maxKNNTaxa= new PluginParameter.Builder<>("maxKNNTaxa",10,Integer.class).range(Range.closed(2,100)).guiName("Preferred haplotype size")
            .description("Preferred haplotype block size in sites (use same as in FILLINFindHaplotypesPlugin)").build();
    private int minLDComparisons=20;

    public LDKNNiImputationPluginEd() {
        super(null, false);
    }

    public LDKNNiImputationPluginEd(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        long time=System.currentTimeMillis();
        System.out.println(time);
        GenotypeTable hetGenotypeTable= ImportUtils.readGuessFormat(hmpFile());
        //GenotypeTable hetGenotypeTable=FilterGenotypeTableBuilder.getInstance(ImportUtils.readGuessFormat(hmpFile())).minorAlleleFreqForSite(0.01,.99).build();
        final GenotypeTable genotypeTable= GenotypeTableBuilder.getHomozygousInstance(hetGenotypeTable);
        System.out.println("genotypeTable.numberOfSites() = " + genotypeTable.numberOfSites());
        System.out.println("genotypeTable.numberOfTaxa() = " + genotypeTable.numberOfTaxa());
        System.out.println(genotypeTable.genotypeAsStringRow(0));


        final int numberOfSites=genotypeTable.numberOfSites();
        final Multimap<Position,Position> highLDMap = ArrayListMultimap.create();
        IntStream.range(0,genotypeTable.numberOfSites()).parallel()
                .forEach(posIndex -> {
                    MinMaxPriorityQueue<LDResult> highestLD=MinMaxPriorityQueue.orderedBy(LDResult.byR2Ordering.reverse())
                            .maximumSize(maxHighLDSites()).create();
                    for (int site2 = 0; site2 < numberOfSites; site2++) {
                        if(posIndex==site2) continue;  //avoid the same site
                        LDResult ld=LinkageDisequilibrium.calculateBitLDForHaplotype(true, minLDComparisons, genotypeTable, posIndex, site2);
                        if(Double.isNaN(ld.r2())) continue;
                        highestLD.add(ld);
                    }
                    List<Position> positionList= new ArrayList<>();
                    for (LDResult result : highestLD) {
                        positionList.add(genotypeTable.positions().get(result.site2()));
                    }
                    highLDMap.putAll(genotypeTable.positions().get(posIndex), positionList);
                });
        System.out.println("LD Map Size:"+highLDMap.keySet().size());
        highLDMap.asMap().entrySet().stream().forEach(es -> System.out.println(es.getKey()+":"+es.getValue().size()));
//        Multiset<Integer> distMapSize=TreeMultiset.create();//highLDMap.keySet().stream().mapToInt(k -> highLDMap.get(k).size()).
//        highLDMap.asMap().entrySet().stream().forEach(en -> distMapSize.add(en.getValue().size()));
//        distMapSize.entrySet().stream()
//                .forEach(System.out::println);

        System.out.println("High LD sites calculated");
        LongAdder genotypesMissing= new LongAdder();
        LongAdder genotypesImputed= new LongAdder();
        GenotypeTableBuilder incSiteBuilder = GenotypeTableBuilder.getSiteIncremental(hetGenotypeTable.taxa());
        IntStream.range(0,hetGenotypeTable.numberOfSites()).parallel().forEach(posIndex -> {
            Position position=hetGenotypeTable.positions().get(posIndex);
            PositionList positionList = PositionListBuilder.getInstance(new ArrayList<>(highLDMap.get(position)));
            byte[] currGenos=hetGenotypeTable.genotypeAllTaxa(posIndex);
            byte[] newGenos=new byte[currGenos.length];
            final GenotypeTable ldGenoTable = GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(genotypeTable, positionList));
            for (int taxon = 0; taxon <currGenos.length; taxon++) {
                if(currGenos[taxon]==GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    genotypesMissing.increment();
                    Multimap<Double, Byte> closeGenotypes=getClosestNonMissingTaxa(hetGenotypeTable.taxa().get(taxon), hetGenotypeTable,
                            ldGenoTable, position, maxKNNTaxa());
                    //System.out.println(closeGenotypes.size());
                    Collection<Byte> genotype=closeGenotypes.get(0.0);  //replace
                    //closeGenotypes.values().size();
                    if(genotype.size()>0) {
                        newGenos[taxon]=genotype.iterator().next();
                        genotypesImputed.increment();
                    } else {
                        newGenos[taxon]=currGenos[taxon];
                    }

                } else {
                    newGenos[taxon]=currGenos[taxon];
                }
            }
            incSiteBuilder.addSite(position,newGenos);
            if(posIndex%10==0) System.out.printf("Position:%d Missing:%d Imputed:%d %n", posIndex, genotypesMissing.intValue(), genotypesImputed.intValue());
            //System.out.println(position.toString()+":"+ibsDistanceMatrix.meanDistance());
        });

        GenotypeTable impGenotypeTable=incSiteBuilder.build();
        ExportUtils.writeToHapmap(impGenotypeTable,outFileBase());

        System.out.println("All Distance matrix calculated");

        return null;
    }

    private Multimap<Double, Byte> getClosestNonMissingTaxa(Taxon inputTaxon, GenotypeTable genotypeTable, GenotypeTable ldGenoTable,
                                                            Position targetPosition, int numberOfTaxa) {
        final int targetPosIdx=genotypeTable.positions().indexOf(targetPosition);
        final int inputTaxonIdx=genotypeTable.taxa().indexOf(inputTaxon);
        byte[] inputTaxonGenotypes=ldGenoTable.genotypeAllSites(inputTaxonIdx);
        MinMaxPriorityQueue<Tuple<Double, Byte>> topTaxa = IntStream.range(0, genotypeTable.numberOfTaxa())
                .filter(closeTaxonIdx -> closeTaxonIdx != inputTaxonIdx)
                .filter(closeTaxonIdx -> genotypeTable.genotype(closeTaxonIdx, targetPosIdx) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)
                .mapToObj(closeTaxonIdx -> new Tuple<>(IBSDistanceMatrix.computeHetDistances(inputTaxonGenotypes, ldGenoTable.genotypeAllSites(closeTaxonIdx), 1)[0], genotypeTable.genotype(closeTaxonIdx, targetPosIdx)))
                .filter(distanceTaxon -> !Double.isNaN(distanceTaxon.x))
                .collect(Collectors.toCollection(() -> MinMaxPriorityQueue.maximumSize(numberOfTaxa).create()));
        if(topTaxa.size()==-1) {
            System.out.println("Crap:"+Arrays.toString(inputTaxonGenotypes));
        }
        final Multimap<Double, Byte> distGenoMap = ArrayListMultimap.create();
        topTaxa.stream().forEach(distGeno -> distGenoMap.put(distGeno.x,distGeno.y));
        return distGenoMap;
    }

    public static double[] distance(byte[] seq1, byte[] seq2, int minSitesCompared) {
        int sameCnt = 0, hetCnt = 0;
        int sites=0;
        for (int i = 0; i < seq1.length; i++) {
            if(seq1[i]==GenotypeTable.UNKNOWN_DIPLOID_ALLELE) continue;
            if(seq2[i]==GenotypeTable.UNKNOWN_DIPLOID_ALLELE) continue;
            sites++;
            if(seq1[i]==seq2[i]) {
                sameCnt++;
            } else if(GenotypeTableUtils.isPartiallyEqual(seq1[i],seq2[i])) {
                hetCnt++;
            }
        }
        double dist=1-(((double)sameCnt+((double)hetCnt/2))/(double)sites);
        if (sites > minSitesCompared) {
            return new double[]{dist, (double) sites};
        } else {
            return new double[]{Double.NaN, (double) sites};
        }

    }

    @Override
    public String getCitation() {
        return "Daniel Money's cool imputation";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "LD based KNN Imputation";
    }

    @Override
    public String getToolTipText() {
        return "Daniel Money's cool imputation";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(LDKNNiImputationPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Input HapMap file of target genotypes to impute. Accepts
     * all file types supported by TASSEL5.
     *
     * @return Target file
     */
    public String hmpFile() {
        return hmpFile.value();
    }

    /**
     * Set Target file. Input HapMap file of target genotypes
     * to impute. Accepts all file types supported by TASSEL5.
     *
     * @param value Target file
     *
     * @return this plugin
     */
    public LDKNNiImputationPluginEd hmpFile(String value) {
        hmpFile = new PluginParameter<>(hmpFile, value);
        return this;
    }

    /**
     * Output file; hmp.txt.gz and .hmp.h5 accepted.
     *
     * @return Output filename
     */
    public String outFileBase() {
        return outFileBase.value();
    }

    /**
     * Set Output filename. Output file; hmp.txt.gz and .hmp.h5
     * accepted.
     *
     * @param value Output filename
     *
     * @return this plugin
     */
    public LDKNNiImputationPluginEd outFileBase(String value) {
        outFileBase = new PluginParameter<>(outFileBase, value);
        return this;
    }

    /**
     * Preferred haplotype block size in sites (use same as
     * in FILLINFindHaplotypesPlugin)
     *
     * @return Maximum number of high LD sites
     */
    public Integer maxHighLDSites() {
        return maxHighLDSites.value();
    }

    /**
     * Set Maximum number of high LD sites. Preferred haplotype
     * block size in sites (use same as in FILLINFindHaplotypesPlugin)
     *
     * @param value Maximum number of high LD sites
     *
     * @return this plugin
     */
    public LDKNNiImputationPluginEd maxHighLDSites(Integer value) {
        maxHighLDSites = new PluginParameter<>(maxHighLDSites, value);
        return this;
    }

    /**
     * Preferred haplotype block size in sites (use same as
     * in FILLINFindHaplotypesPlugin)
     *
     * @return Preferred haplotype size
     */
    public Integer maxKNNTaxa() {
        return maxKNNTaxa.value();
    }

    /**
     * Set Preferred haplotype size. Preferred haplotype block
     * size in sites (use same as in FILLINFindHaplotypesPlugin)
     *
     * @param value Preferred haplotype size
     *
     * @return this plugin
     */
    public LDKNNiImputationPluginEd maxKNNTaxa(Integer value) {
        maxKNNTaxa = new PluginParameter<>(maxKNNTaxa, value);
        return this;
    }
}
