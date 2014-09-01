package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.analysis.gbs.ProcessFastQFile;
import net.maizegenetics.analysis.gbs.TagDistributionMap;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.util.DirectoryCrawler;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Develops a discovery TBT file from a set of GBS sequence files.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 *
 */
public class DiscoveryTBTPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(DiscoveryTBTPlugin.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing FASTQ files in text or gzipped text.\n"
                    + "     NOTE: Directory will be searched recursively and should\n"
                    + "     be written WITHOUT a slash after its name.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myEnzyme = new PluginParameter.Builder<>("e", null, String.class).guiName("Enzyme").required(true)
            .description("Enzyme used to create the GBS library, if it differs from the one listed in the key file").build();
    private PluginParameter<Integer> myMinTagCount = new PluginParameter.Builder<>("c", 10, Integer.class).guiName("Min Tag Count")
            .description("Minimum tag count").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<>("o", null, String.class).guiName("Output TBT File").required(true).outFile()
            .description("Output TBT file").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();
    private PluginParameter<Integer> myMaxMapMemoryInMb = new PluginParameter.Builder<>("mxMapMem", 8000, Integer.class).guiName("Maximum Map Memory in Mb").required(false)
            .description("Maximum size for the tag distribution map in Mb").build();

    private int maxMapSize=20_000_000;
    private TagDistributionMap tagCntMap;

    private static final String inputFileGlob="glob:*{.fq,fq.gz,fastq,fastq.txt,fastq.gz,fastq.txt.gz,_sequence.txt,_sequence.txt.gz}";
    private static final int keyFileTaxaNameIndex=3;
    private static final int keyFilePrepIDIndex=7;
    static final String sampleNameField="FullSampleName";




    public DiscoveryTBTPlugin() {
        super(null, false);
    }

    public DiscoveryTBTPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        tagCntMap=new TagDistributionMap(maxMapSize,(long)myMaxMapMemoryInMb.value()*1_000_000l);
        try {
            //Get the map of taxa indices from the key files
            //Get the list of fastq files
            Path keyPath= Paths.get(keyFile()).toAbsolutePath();
            java.util.List<Path> inputSeqFiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
            if(inputSeqFiles.isEmpty()) {
                myLogger.warn("No files matching:"+inputFileGlob);
                return null;
            }
            TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), sampleNameField, new HashMap<String, String>(), true);
            //setup
            int numThreads=Runtime.getRuntime().availableProcessors();
            ExecutorService pool= Executors.newFixedThreadPool(numThreads - 1);
            for (Path inputSeqFile : inputSeqFiles) {
                System.out.println("Using File:"+inputSeqFile.toString());
                ProcessFastQFile pb=new ProcessFastQFile(masterTaxaList,keyPath, inputSeqFile, enzyme(),
                        minimumQualityScore(), tagCntMap);
                //pb.run();
                pool.execute(pb);
                //countTags(myKeyFile.value(), myEnzyme.value(), inputSeqFile);
                System.out.println("tagCntMap.size()="+tagCntMap.size());
            }
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BuilderFromHapMap: processing threads timed out.");
            }
            System.out.println("tagCntMap.size()="+tagCntMap.size());
            System.out.println("Memory Size"+tagCntMap.estimateMapMemorySize());
            tagCntMap.removeTagByCount(0);
            System.out.println("tagCntMap.size()="+tagCntMap.size());
            System.out.println("Memory Size"+tagCntMap.estimateMapMemorySize());
            int cnt=0;
            int totalmatch=0;
            for (Map.Entry<Tag, TaxaDistribution> entry : tagCntMap.entrySet()) {
                if (cnt<50) System.out.println(entry.getKey().toString()+":"+entry.getValue().toString());
                cnt++;
                totalmatch+=entry.getValue().totalDepth();
            }
            System.out.println("totalmatch = " + totalmatch);
            TagDataWriter tdw=new TagDataSQLite(myOutputFile.value());
            tdw.putTaxaList(masterTaxaList);
            tdw.putAllTag(tagCntMap.keySet());
            tdw.putTaxaDistribution(tagCntMap);
            //TagsByTaxaHDF5Builder.create(myOutputFile.value(), tagCntMap, masterTaxaList);

        } catch(Exception e) {
            e.printStackTrace();
        }

        return new DataSet(new Datum("TagMap",tagCntMap,""),this);
    }



    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    public static void main(String[] args) {
         GeneratePluginCode.generate(DiscoveryTBTPlugin.class);
    }

    /**
     * Input directory containing FASTQ files in text or gzipped
     * text.
     *      NOTE: Directory will be searched recursively and
     * should
     *      be written WITHOUT a slash after its name.
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDir.value();
    }

    /**
     * Set Input Directory. Input directory containing FASTQ
     * files in text or gzipped text.
     *      NOTE: Directory will be searched recursively and
     * should
     *      be written WITHOUT a slash after its name.
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin inputDirectory(String value) {
        myInputDir = new PluginParameter<>(myInputDir, value);
        return this;
    }

    /**
     * Key file listing barcodes distinguishing the samples
     *
     * @return Key File
     */
    public String keyFile() {
        return myKeyFile.value();
    }

    /**
     * Set Key File. Key file listing barcodes distinguishing
     * the samples
     *
     * @param value Key File
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }

    /**
     * Enzyme used to create the GBS library, if it differs
     * from the one listed in the key file
     *
     * @return Enzyme
     */
    public String enzyme() {
        return myEnzyme.value();
    }

    /**
     * Set Enzyme. Enzyme used to create the GBS library,
     * if it differs from the one listed in the key file
     *
     * @param value Enzyme
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

    /**
     * Minimum tag count
     *
     * @return Min Tag Count
     */
    public Integer minTagCount() {
        return myMinTagCount.value();
    }

    /**
     * Set Min Tag Count. Minimum tag count
     *
     * @param value Min Tag Count
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin minTagCount(Integer value) {
        myMinTagCount = new PluginParameter<>(myMinTagCount, value);
        return this;
    }

    /**
     * Output directory to contain .cnt files (one per FASTQ
     * file (one per FASTQ file)
     *
     * @return Output Directory
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output Directory. Output directory to contain .cnt
     * files (one per FASTQ file (one per FASTQ file)
     *
     * @param value Output Directory
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<String>(myOutputFile, value);
        return this;
    }

    /**
     * Minimum quality score within the barcode and read length
     * to be accepted
     *
     * @return Minimum quality score
     */
    public Integer minimumQualityScore() {
        return myMinQualScore.value();
    }

    /**
     * Set Minimum quality score. Minimum quality score within
     * the barcode and read length to be accepted
     *
     * @param value Minimum quality score
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }

    /**
     * Maximum size for the tag distribution map in Mb
     * @return Maximum Map Memory in Mb
     */
    public Integer maximumMapMemoryInMb() {
        return myMaxMapMemoryInMb.value();
    }

    /**
     * Set Maximum Map Memory in Mb. Maximum size for the tag distribution map in Mb
     * @param value Maximum Map Memory in Mb
     * @return this plugin
     */
    public DiscoveryTBTPlugin maximumMapMemoryInMb(Integer value) {
        myMaxMapMemoryInMb = new PluginParameter<>(myMaxMapMemoryInMb, value);
        return this;
    }


    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Discovery Tags By Taxa";
    }

    @Override
    public String getToolTipText() {
        return "Discovery Tags By Taxa";
    }
}
