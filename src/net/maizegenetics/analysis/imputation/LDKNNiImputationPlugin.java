package net.maizegenetics.analysis.imputation;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.MinMaxPriorityQueue;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import java.util.stream.IntStream;
import javax.swing.ImageIcon;

import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;

/**
 * Need to fill in this
 *
 * @author Daniel Money
 * @author Ed Buckler
 */
public class LDKNNiImputationPlugin extends AbstractPlugin {

    private PluginParameter<String> hmpFile = new PluginParameter.Builder<>("i", null, String.class)
            .guiName("Target file")
            .inFile()
            .required(true)
            .description("Input HapMap file of target genotypes to impute. Accepts all file types supported by TASSEL5.")
            .build();

    private PluginParameter<String> outFileBase = new PluginParameter.Builder<>("o", null, String.class)
            .guiName("Output filename")
            .outFile()
            .required(true)
            .description("Output file; hmp.txt.gz and .hmp.h5 accepted.")
            .build();

    private PluginParameter<Integer> maxHighLDSites = new PluginParameter.Builder<>("maxHighLDSSites", 20, Integer.class)
            .range(Range.closed(2, 2000))
            .guiName("Max High LD Sites")
            .description("Maximum number of sites in high LD to use in imputation")
            .build();

    private PluginParameter<Integer> maxKNNTaxa = new PluginParameter.Builder<>("maxKNNTaxa", 20, Integer.class)
            .range(Range.closed(2, 200))
            .guiName("Max kNN taxa")
            .description("Maximum number of neighbours to use in imputation")
            .build();

    public LDKNNiImputationPlugin() {
        super(null, false);
    }

    public LDKNNiImputationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        // Load in the genotype table
        GenotypeTable hetGenotypeTable = ImportUtils.readGuessFormat(hmpFile());
        // Create a copy of the table with hets removed for use in LD calcuations
        final GenotypeTable genotypeTable = GenotypeTableBuilder.getHomozygousInstance(hetGenotypeTable);

        // Create a multimap of the snps in higest LD with each SNP
        Multimap<Position, Position> highLDMap = getHighLDMap(genotypeTable, maxHighLDSites());

        LongAdder genotypesMissing = new LongAdder();
        LongAdder genotypesImputed = new LongAdder();
        GenotypeTableBuilder incSiteBuilder = GenotypeTableBuilder.getSiteIncremental(hetGenotypeTable.taxa());
        IntStream.range(0, hetGenotypeTable.numberOfSites()).parallel().forEach(posIndex ->
        {
            Position position = hetGenotypeTable.positions().get(posIndex);
            PositionList positionList = PositionListBuilder.getInstance(new ArrayList<>(highLDMap.get(position)));
            byte[] currGenos = hetGenotypeTable.genotypeAllTaxa(posIndex);
            byte[] newGenos = new byte[currGenos.length];
            final GenotypeTable ldGenoTable = GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(genotypeTable, positionList));
            for (int taxon = 0; taxon < currGenos.length; taxon++) {
                if (currGenos[taxon] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    genotypesMissing.increment();

                    if (isInvariant(currGenos)) {
                        newGenos[taxon] = getInvariantGeno(currGenos);
                        genotypesImputed.increment();
                    } else {

                        Multimap<Double, Byte> closeGenotypes = getClosestNonMissingTaxa(hetGenotypeTable.taxa().get(taxon), hetGenotypeTable,
                                ldGenoTable, position, maxKNNTaxa());
                        if (closeGenotypes.isEmpty())
                        //Collection<Byte> genotype=closeGenotypes.get(0.0);  //replace
                        //if(genotype.size()==0)
                        {
                            newGenos[taxon] = currGenos[taxon];
                        } else {
                            newGenos[taxon] = impute(closeGenotypes, 20);
                            genotypesImputed.increment();
                        }
                    }
                } else {
                    newGenos[taxon] = currGenos[taxon];
                }

            }

            incSiteBuilder.addSite(position, newGenos);
            if (posIndex % 10 == 0) {
                System.out.printf("Position:%d Missing:%d Imputed:%d %n", posIndex, genotypesMissing.intValue(), genotypesImputed.intValue());
            }
            //System.out.println(position.toString()+":"+ibsDistanceMatrix.meanDistance());
        });
        GenotypeTable impGenotypeTable = incSiteBuilder.build();
        ExportUtils.writeToHapmap(impGenotypeTable, outFileBase());

        System.out.println("All Distance matrix calculated");

        return null;
    }

    private boolean isInvariant(byte[] genos) {
        byte foundGeno = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;

        for (byte geno : genos) {
            if (geno != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                if (foundGeno == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    foundGeno = geno;
                } else {
                    if (foundGeno != geno) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    private byte getInvariantGeno(byte[] genos) {
        for (byte geno : genos) {
            if (geno != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                return geno;
            }
        }

        return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
    }

    private Multimap<Double, Byte> getClosestNonMissingTaxa(Taxon inputTaxon, GenotypeTable genotypeTable, GenotypeTable ldGenoTable,
                                                            Position targetPosition, int numberOfTaxa) {
        final int targetPosIdx = genotypeTable.positions().indexOf(targetPosition);
        final int inputTaxonIdx = genotypeTable.taxa().indexOf(inputTaxon);
        byte[] inputTaxonGenotypes = ldGenoTable.genotypeAllSites(inputTaxonIdx);
        MinMaxPriorityQueue<Tuple<Double, Byte>> topTaxa = IntStream.range(0, genotypeTable.numberOfTaxa())
                .filter(closeTaxonIdx -> closeTaxonIdx != inputTaxonIdx)
                .filter(closeTaxonIdx -> genotypeTable.genotype(closeTaxonIdx, targetPosIdx) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)
                .mapToObj(closeTaxonIdx -> new Tuple<>(IBSDistanceMatrix.computeHetDistances(inputTaxonGenotypes, ldGenoTable.genotypeAllSites(closeTaxonIdx), 10)[0], genotypeTable.genotype(closeTaxonIdx, targetPosIdx)))
                .filter(distanceTaxon -> !Double.isNaN(distanceTaxon.x))
                .collect(Collectors.toCollection(() -> MinMaxPriorityQueue.maximumSize(numberOfTaxa).create()));
        final Multimap<Double, Byte> distGenoMap = ArrayListMultimap.create();
        topTaxa.stream().forEach(distGeno -> distGenoMap.put(distGeno.x, distGeno.y));
        return distGenoMap;
    }

    private Multimap<Position, Position> getHighLDMap(GenotypeTable genotypeTable, int numberOfSNPs) {
        Multimap<Position, Position> highLDMap = ArrayListMultimap.create();

        final int numberOfSites = genotypeTable.numberOfSites();

        IntStream.range(0, genotypeTable.numberOfSites()).parallel()
                .forEach(posIndex ->
                {
                    MinMaxPriorityQueue<LDResult> highestLD = MinMaxPriorityQueue.orderedBy(LDResult.byR2Ordering.reverse())
                            .maximumSize(numberOfSNPs).create();
                    for (int site2 = 0; site2 < numberOfSites; site2++) {
                        if (posIndex == site2) {
                            continue;
                        }
                        LDResult ld = LinkageDisequilibrium.calculateBitLDForHaplotype(true, 20, genotypeTable, posIndex, site2);
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
                });

        return highLDMap;
    }

    private byte impute(Multimap<Double, Byte> distGeno, int useLDSites) {
        // Create an array to store the weighted counts
        double[] weightedCount = new double[255];

        // For each distance to genotype / genotype pair update the weighted counts
        distGeno.entries().forEach(entry ->
        {
            weightedCount[entry.getValue() + 128] += 1.0 / (1.0 + useLDSites * entry.getKey());
        });

        // Find the best genotype
        int bestGeno = 0;
        double bestWeightedCount = weightedCount[0];
        for (int i = 1; i < 255; i++) {
            if (weightedCount[i] > bestWeightedCount) {
                bestWeightedCount = weightedCount[i];
                bestGeno = i;
            }
        }

        return (byte) (bestGeno - 128);
    }

    @Override
    public String getCitation() {
        return "";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }

    @Override
    public String getToolTipText() {
        return null;
    }

    public static void main(String[] args) {
        GeneratePluginCode.generate(LDKNNiImputationPlugin.class);
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
    // TODO: Replace <Type> with specific type.
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Input HapMap file of target genotypes to impute. Accepts all file types
     * supported by TASSEL5.
     *
     * @return Target file
     */
    public String hmpFile() {
        return hmpFile.value();
    }

    /**
     * Set Target file. Input HapMap file of target genotypes to impute. Accepts
     * all file types supported by TASSEL5.
     *
     * @param value Target file
     * @return this plugin
     */
    public LDKNNiImputationPlugin hmpFile(String value) {
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
     * Set Output filename. Output file; hmp.txt.gz and .hmp.h5 accepted.
     *
     * @param value Output filename
     * @return this plugin
     */
    public LDKNNiImputationPlugin outFileBase(String value) {
        outFileBase = new PluginParameter<>(outFileBase, value);
        return this;
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
     * Set Max High LD Sites. Maximum number of sites in high LD to use in
     * imputation
     *
     * @param value Max High LD Sites
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
     * Set Max kNN taxa. Maximum number of neighbours to use in imputation
     *
     * @param value Max kNN taxa
     * @return this plugin
     */
    public LDKNNiImputationPlugin maxKNNTaxa(Integer value) {
        maxKNNTaxa = new PluginParameter<>(maxKNNTaxa, value);
        return this;
    }
}
