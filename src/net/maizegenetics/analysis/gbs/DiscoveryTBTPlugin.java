package net.maizegenetics.analysis.gbs;

import com.google.common.collect.BiMap;
import com.google.common.collect.ImmutableMap;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.nio.file.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

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
            Path keyPath=Paths.get(keyFile()).toAbsolutePath();
            List<Path> inputSeqFiles=DirectoryCrawler.listPaths(inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
            if(inputSeqFiles.isEmpty()) {
                myLogger.warn("No files matching:"+inputFileGlob);
                return null;
            }
            TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(),sampleNameField,new HashMap<String, String>(),true);
            //setup
            int numThreads=Runtime.getRuntime().availableProcessors();
            ExecutorService pool= Executors.newFixedThreadPool(numThreads-1);
            for (Path inputSeqFile : inputSeqFiles) {
                System.out.println("Using File:"+inputSeqFile.toString());
                ProcessFastQFile pb=new ProcessFastQFile(masterTaxaList,keyPath, inputSeqFile, enzyme(),
                        minimumQualityScore(), tagCntMap);
                //pb.run();
                pool.execute(pb);
                //countTags(myKeyFile.value(), myEnzyme.value(), inputSeqFile);\
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
            TagsByTaxaHDF5Builder.create(myOutputFile.value(), tagCntMap, masterTaxaList);

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
        myOutputFile = new PluginParameter<>(myOutputFile, value);
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


/**
 * This concurrnet HashMap constrain the size of the map, and purges low distribution count tags when the size needs
 * to be reduced.
 *
 *
 */
class TagDistributionMap extends ConcurrentHashMap<Tag,TaxaDistribution> {
    private final long maxMemorySize;
    private int minDepthToRetainInMap=1;
    private AtomicLong putCntSinceMemoryCheck=new AtomicLong(0L);  //since we don't need an exact count of the puts,
    //this may be overkill
    private long checkFreq;  //numbers of puts to check size

    TagDistributionMap(int initialCapacity, long maxMemorySize) {
        super(initialCapacity);
        this.maxMemorySize=maxMemorySize;
        checkFreq=(maxMemorySize/(100L*10L));  //this is checking roughly to keep the map within 10% of the max
    }

    @Override
    public TaxaDistribution put(Tag key, TaxaDistribution value) {
        if(putCntSinceMemoryCheck.incrementAndGet()>checkFreq) {
            checkMemoryAndReduceIfNeeded();
            putCntSinceMemoryCheck.set(0);
        }
        return super.put(key, value);
    }

    public synchronized void checkMemoryAndReduceIfNeeded() {
        System.out.println("TagDistributionMap.checkMemoryAndReduceIfNeeded"+estimateMapMemorySize());
        while(estimateMapMemorySize()>maxMemorySize) {
            reduceMapSize();
        }
    }


    public synchronized void removeTagByCount(int minCnt) {
        for (Tag tag : keySet()) {
            if(get(tag).totalDepth()<minCnt) remove(tag);
        }
    }

    private synchronized void reduceMapSize() {
        System.out.println("reduceMapSize()"+minDepthToRetainInMap);
        removeTagByCount(minDepthToRetainInMap);
        minDepthToRetainInMap++;
    }

    public long estimateMapMemorySize() {
        long size=0;
        int cnt=0;
        for (Map.Entry<Tag, TaxaDistribution> entry : entrySet()) {
            size+=8+16+1; //Tag size
            size+=16; //Map references
            size+=entry.getValue().memorySize();
            cnt++;
            if(cnt>10000) break;
        }
        long estSize=(size()/cnt)*size;
        return estSize;
    }
}


class ProcessFastQFile implements Runnable {
    private static final Logger myLogger = Logger.getLogger(DiscoveryTBTPlugin.class);
    private final TaxaList masterTaxaList;
    private final Path keyPath;
    private final Path fastQPath;
    private final Map<Tag,TaxaDistribution> masterTagTaxaMap;
    private String enzyme;
    private final int minQuality;
    private final String sampleNameField=DiscoveryTBTPlugin.sampleNameField;
    private final String enzymeField="Enzyme";
    private final String flowcellField="Flowcell";
    private final String laneField="Lane";



    ProcessFastQFile(TaxaList masterTaxaList, Path keyPath, Path fastQPath, String enzyme,
                             int minQuality, Map<Tag,TaxaDistribution> masterTagTaxaMap) {
        this.masterTaxaList=masterTaxaList;
        this.keyPath = keyPath;
        this.fastQPath =fastQPath;
        this.masterTagTaxaMap=masterTagTaxaMap;
        this.enzyme=enzyme;
        this.minQuality=minQuality;
    }

    @Override
    public void run() {
        TaxaList tl=getLaneAnnotatedTaxaList(keyPath, fastQPath);
        if(enzyme==null) enzyme=tl.get(0).getTextAnnotation(enzymeField)[0];
        ParseBarcodeRead2 thePBR=new ParseBarcodeRead2(tl, enzyme, masterTaxaList);
        processFastQ(fastQPath,thePBR);
    }

    private TaxaList getLaneAnnotatedTaxaList(Path keyPath, Path fastQpath) {
        //myLogger.info("Reading FASTQ file: " + fastqFile.toString());
        String[] filenameField = fastQpath.getFileName().toString().split("_");
        TaxaList annoTL;
        ParseBarcodeRead2 thePBR;  // this reads the key file and store the expected barcodes for this lane
        if (filenameField.length == 3) {
            annoTL = TaxaListIOUtils.readTaxaAnnotationFile(keyPath.toAbsolutePath().toString(), sampleNameField,
                    ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[1]), false);
        } else if (filenameField.length == 4) {
            annoTL = TaxaListIOUtils.readTaxaAnnotationFile(keyPath.toAbsolutePath().toString(),sampleNameField,
                    ImmutableMap.of(flowcellField, filenameField[0], laneField, filenameField[2]),false);
        }
        else if (filenameField.length == 5) {
            annoTL = TaxaListIOUtils.readTaxaAnnotationFile(keyPath.toAbsolutePath().toString(),sampleNameField,
                    ImmutableMap.of(flowcellField, filenameField[1], laneField, filenameField[3]),false);
        } else {
            myLogger.error("Error in parsing file name: " + fastQpath.toString());
            myLogger.error("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
            myLogger.error("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
            return null;
        }
        return annoTL;
    }

    private void processFastQ(Path fastqFile, ParseBarcodeRead2 thePBR) {
            int goodBarcodedReads = 0;
            int maxTaxaNumber=masterTaxaList.size();

            try {
                BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1<<22);
                int currLine = 0;
                long time=System.nanoTime();
                int allReads = 0;
                goodBarcodedReads = 0;
                String sequence = "";
                String qualityScore = "";
                String temp = br.readLine();
                while (temp != null) {
                    currLine++;
                    if(currLine%(1<<20)==0) {
                        System.out.printf("File line %d processing rate %d ns/line  %n",currLine, (System.nanoTime()-time)/currLine);
                    }
                    try {
                        //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
                        if ((currLine + 2) % 4 == 0) {
                            sequence = temp;
                        } else if (currLine % 4 == 0) {
                            //todo look at qualities and only process good ones
                            qualityScore = temp;
                            allReads++;
                            //After quality score is read, decode barcode using the current sequence & quality  score
                            ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sequence, qualityScore, true, minQuality);
                            if (rr != null) {
                                goodBarcodedReads++;
                                Tag tg=TagBuilder.instance(rr.getRead(), rr.getLength());
                                TaxaDistribution iC=masterTagTaxaMap.get(tg);
                                if(iC==null) {
                                    masterTagTaxaMap.put(tg,TaxaDistBuilder.create(maxTaxaNumber,rr.getTaxonIndex()));
                                }
                                else if(iC.totalDepth()==1) {
                                    masterTagTaxaMap.put(tg,TaxaDistBuilder.create(iC).increment(rr.getTaxonIndex()));
                                } else {
                                   iC.increment(rr.getTaxonIndex());
                                }
                            }
                            if (allReads % 1000000 == 0) {
                                myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads);
                            }
                        }
                    } catch (NullPointerException e) {
                        e.printStackTrace();
                        myLogger.error("Unable to correctly parse the sequence and: " + sequence
                                + " and quality score: " + qualityScore + " from fastq file.  Your fastq file may have been corrupted.");
                        System.exit(1);
                    }
                    temp = br.readLine();
                }

                myLogger.info("Total number of reads in lane=" + allReads);
                myLogger.info("Total number of good barcoded reads=" + goodBarcodedReads);
                myLogger.info("Timing process (sorting, collapsing, and writing TagCount to file).");
                long timePoint1 = System.currentTimeMillis();
                System.out.println("tagCntMap size"+masterTagTaxaMap.size());
//                for (Map.Entry<Tag, TaxaDistribution> entry : masterTagTaxaMap.entrySet()) {
//                    if(entry.getKey().seqLength()==64 && entry.getValue().totalDepth()>1000) {
//                        System.out.println(entry.toString());
//                        System.out.println(entry.getValue().taxaDepthMap());
//                    }
//
//                }

                //theTC.writeTagCountFile(outputDir + File.separator + countFileNames[laneNum], FilePacking.Byte, minCount);
                myLogger.info("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
                br.close();


            } catch (Exception e) {
                myLogger.error("Good Barcodes Read: " + goodBarcodedReads);
                e.printStackTrace();
            }
    }

}
