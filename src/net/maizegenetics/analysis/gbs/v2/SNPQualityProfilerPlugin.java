package net.maizegenetics.analysis.gbs.v2;


import com.google.common.collect.Multimap;
import com.google.common.primitives.Ints;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * Created by edbuckler on 11/21/14.
 */
public class SNPQualityProfilerPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(SAMToGBSdbPlugin.class);

    private PluginParameter<String> myTaxaFile = new PluginParameter.Builder<String>("taxa", null, String.class).guiName("Taxa List File").required(true).inFile()
            .description("Name of taxa list input file in taxa list format").build();
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("GBS DB File").required(true).outFile()
            .description("Name of output file (e.g. GBSv2.db)").build();
    private TagDataSQLite tagDataWriter;

    public SNPQualityProfilerPlugin() {
        super(null, false);
    }

    public SNPQualityProfilerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        tagDataWriter =new TagDataSQLite(gBSDBFile());
        TaxaList taxaList=tagDataWriter.getTaxaList();
        TaxaList subTaxa=TaxaListIOUtils.readTaxaAnnotationFile(taxaListFile(),"<NAME>");
        System.out.println("sublist");


        subTaxa.stream().filter(t -> taxaList.indexOf(t)<0).forEach(t-> System.err.println("Missing taxon from master:" + t));
        int[] subsetIndices=subTaxa.stream().mapToInt(taxaList::indexOf).filter(i -> i > -1).sorted().toArray();
       // int[] subsetIndices= IntStream.range(0,taxaList.numberOfTaxa()).toArray();  //for testing using all taxa

        tagDataWriter.getSNPPositions().stream()
                .forEach(position -> {
                    Multimap<Allele,TaxaDistribution> aTDMMap=tagDataWriter.getAllelesTaxaDistForSNP(position);
                    Map<Allele, int[]> subDepths=convertToSubsetMap(aTDMMap,subsetIndices);
                    //filter only for Major and minor

                    calcMinorAlleleDepth(subDepths);
                    calcMAF(subDepths);
                    System.out.println(subDepths.toString());
                    for (Map.Entry<Allele, int[]> alleleEntry : subDepths.entrySet()) {
                        //System.out.println(alleleEntry.getKey().toString()+":"+Arrays.toString(alleleEntry.getValue()));
                        //System.out.println(alleleEntry.getKey().toString()+":"+ Arrays.stream(alleleEntry.getValue()).sum());
                    }
                });


        System.out.println(Arrays.toString(subsetIndices));



        //create subset array to make fast calculations
        //grab a stream of sites and taxa distributions
        //create subset vector that make taxa used
        //process in parallel


        //TaxaListIOUtils.exportAnnotatedTaxaListTable();
        return null;
    }

    /*
    Reduces the TaxaDistributions to depths by allele-> selected subset of taxa
    Some alleles are represented by multiple tags and they are collapsed.
     */
    private Map<Allele, int[]> convertToSubsetMap(Multimap<Allele,TaxaDistribution> aTDMMap, int[] subsetIndices) {
        Map<Allele, int[]> result=new HashMap<>();
        for (Allele allele : aTDMMap.keySet()) {
            int[] subDepths=new int[subsetIndices.length];
            for (TaxaDistribution taxaDistribution : aTDMMap.get(allele)) {
                System.out.println("TD:"+allele.toString()+":"+taxaDistribution.totalDepth());
                int[] depths=taxaDistribution.depths();
                for (int i = 0; i < subDepths.length; i++) {
                    subDepths[i]+=depths[subsetIndices[i]];
                }
            }
            result.put(allele,subDepths);
        }
        return result;
    }

    private double calcMAF(Map<Allele, int[]> alleleDepthMap) {
        int[] sum=alleleDepthMap.values().stream()
                .mapToInt(depths -> (int)Arrays.stream(depths).filter(d -> d>0).count())
                .sorted().toArray();
        System.out.println("MAC:"+Arrays.toString(sum));
        return (sum.length)>1?sum[sum.length-2]:0;
    }

    //Minor allele depth
    private double calcMinorAlleleDepth(Map<Allele, int[]> alleleDepthMap) {
        int[] sum=alleleDepthMap.values().stream()
                    .mapToInt(depths -> Arrays.stream(depths).sum())
                    .sorted().toArray();
        System.out.println("MAD:"+Arrays.toString(sum));
        return (sum.length)>1?sum[sum.length-2]:0;
    }


//    The following getters and setters were auto-generated.
//    Please use this method to re-generate.
//
//    public static void main(String[] args) {
//         GeneratePluginCode.generate(SNPQualityProfilerPlugin.class);
//    }


    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "SNP Quality Profiler";
    }

    @Override
    public String getToolTipText() {
        return "SNP Quality Profiler";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(SNPQualityProfilerPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TagData runPlugin(DataSet input) {
        return (TagData) performFunction(input).getData(0).getData();
    }

    /**
     * Name of taxa list input file in taxa list format
     *
     * @return Taxa List File
     */
    public String taxaListFile() {
        return myTaxaFile.value();
    }

    /**
     * Set Taxa List File. Name of taxa list input file in
     * taxa list format
     *
     * @param value Taxa List File
     *
     * @return this plugin
     */
    public SNPQualityProfilerPlugin taxaListFile(String value) {
        myTaxaFile = new PluginParameter<>(myTaxaFile, value);
        return this;
    }

    /**
     * Name of output file (e.g. GBSv2.db)
     *
     * @return GBS DB File
     */
    public String gBSDBFile() {
        return myDBFile.value();
    }

    /**
     * Set GBS DB File. Name of output file (e.g. GBSv2.db)
     *
     * @param value GBS DB File
     *
     * @return this plugin
     */
    public SNPQualityProfilerPlugin gBSDBFile(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }
}
