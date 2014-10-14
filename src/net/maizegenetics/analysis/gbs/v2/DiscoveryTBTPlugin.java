package net.maizegenetics.analysis.gbs.v2;

import com.google.common.collect.ImmutableMap;
import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;


import javax.swing.ImageIcon;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;

/**
 * Develops a discovery TBT file from a set of GBS sequence files.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 *
 * @author Ed Buckler
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
    private PluginParameter<Integer> myMaxTagLength = new PluginParameter.Builder<>("mxTagL", 64, Integer.class).guiName("Maximum Tag Length")
            .description("Maximum Tag Length").build();
    private PluginParameter<Integer> myMinTagLength = new PluginParameter.Builder<>("mnTagL", 20, Integer.class).guiName("Minimum Tag Length")
            .description("Minimum Tag Length").build();
    private PluginParameter<Integer> myMinTagCount = new PluginParameter.Builder<>("c", 10, Integer.class).guiName("Min Tag Count")
            .description("Minimum tag count").build();
    private PluginParameter<String> myOutputDB = new PluginParameter.Builder<>("o", null, String.class).guiName("Output Database File").required(true).outFile()
            .description("Output Database File").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();
    private PluginParameter<Integer> myMaxMapMemoryInMb = new PluginParameter.Builder<>("mxMapMem", 8000, Integer.class).guiName("Maximum Map Memory in Mb").required(false)
            .description("Maximum size for the tag distribution map in Mb").build();

    private TagDistributionMap tagCntMap;

    static final String inputFileGlob="glob:*{.fq,fq.gz,fastq,fastq.txt,fastq.gz,fastq.txt.gz,_sequence.txt,_sequence.txt.gz}";
    static final String sampleNameField="FullSampleName";
    static final String flowcellField="Flowcell";
    static final String laneField="Lane";
    static final String barcodeField="Barcode";

    public DiscoveryTBTPlugin() {
        super(null, false);
    }

    public DiscoveryTBTPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        tagCntMap=new TagDistributionMap(myMaxMapMemoryInMb.value()/100,(long)myMaxMapMemoryInMb.value()*1_000_000l);
        try {
            //Get the list of fastq files
            Path keyPath= Paths.get(keyFile()).toAbsolutePath();
            List<Path> inputSeqFiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
            if(inputSeqFiles.isEmpty()) {
                myLogger.warn("No files matching:"+inputFileGlob);
                return null;
            }
            TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), sampleNameField, new HashMap<>(), true);
            inputSeqFiles.parallelStream()
                    .forEach(inputSeqFile -> {
                        processFastQFile(masterTaxaList,keyPath, inputSeqFile, enzyme(),
                                minimumQualityScore(), tagCntMap, maximumTagLength());
                    });
            removeSecondCutSitesFromMap(new GBSEnzyme(enzyme()));
            tagCntMap.removeTagByCount(myMinTagCount.value());

            TagDataWriter tdw=new TagDataSQLite(myOutputDB.value());
            tdw.putTaxaList(masterTaxaList);
            tdw.putAllTag(tagCntMap.keySet());
            tdw.putTaxaDistribution(tagCntMap);
            ((TagDataSQLite)tdw).close();  //todo autocloseable should do this but it is not working.

        } catch(Exception e) {
            e.printStackTrace();
        }

        return new DataSet(new Datum("TagMap",tagCntMap,""),this);
    }

    private void processFastQFile(TaxaList masterTaxaList, Path keyPath, Path fastQPath, String enzymeName,
                     int minQuality, Map<Tag,TaxaDistribution> masterTagTaxaMap, int preferredTagLength) {
        TaxaList tl=getLaneAnnotatedTaxaList(keyPath, fastQPath);
        BarcodeTrie barcodeTrie=initializeBarcodeTrie(tl, masterTaxaList, new GBSEnzyme(enzymeName));
        processFastQ(fastQPath,barcodeTrie,masterTaxaList,masterTagTaxaMap,preferredTagLength,minQuality);
    }

    /**
     * Produces a trie for sorting the read
     * @param taxaList the taxaList of the current flowcell lanes that is annotated with barcode information
     * @param masterTaxaList  the mastertaxaList provides the taxaIndex
     * @param myEnzyme
     * @return Barcode trie for examining the prefixes
     */
    static BarcodeTrie initializeBarcodeTrie(TaxaList taxaList, TaxaList masterTaxaList, GBSEnzyme myEnzyme){
        BarcodeTrie aTrie=new BarcodeTrie();
        for (Taxon taxon : taxaList) {
            int masterIndex=masterTaxaList.indexOf(taxon.getName());
            Barcode theBC = new Barcode(taxon.getTextAnnotation(barcodeField)[0], myEnzyme.initialCutSiteRemnant(), taxon.getName(),
                    masterIndex,taxon.getTextAnnotation(flowcellField)[0],taxon.getTextAnnotation("Lane")[0]);
            aTrie.addBarcode(theBC);
        }
        return aTrie;
    }

    /**
     * Returns an annotated taxaList based on a Keyfile for GBS
     * @param keyPath
     * @param fastQpath
     * @return
     */
    static TaxaList getLaneAnnotatedTaxaList(Path keyPath, Path fastQpath) {
        String[] filenameField = fastQpath.getFileName().toString().split("_");
        TaxaList annoTL;
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

    private void processFastQ(Path fastqFile, BarcodeTrie barcodeTrie, TaxaList masterTaxaList,
                              Map<Tag,TaxaDistribution> masterTagTaxaMap, int preferredTagLength, int minQual) {
        int allReads=0, goodBarcodedReads = 0;
        int maxTaxaNumber=masterTaxaList.size();
        try {
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;
            while ((seqAndQual=readFastQBlock(br,allReads)) != null) {
                allReads++;
                //After quality score is read, decode barcode using the current sequence & quality  score
                Barcode barcode=barcodeTrie.longestPrefix(seqAndQual[0]);
                if(barcode==null) continue;
                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual)<(barcode.getBarLength()+preferredTagLength)) continue;
                }
                Tag tag=TagBuilder.instance(seqAndQual[0].substring(barcode.getBarLength(), barcode.getBarLength()+preferredTagLength)).build();
                if(tag==null) continue;   //null occurs when any base was not A, C, G, T
                goodBarcodedReads++;
                TaxaDistribution taxaDistribution=masterTagTaxaMap.get(tag);
                if(taxaDistribution==null) {
                    masterTagTaxaMap.put(tag,TaxaDistBuilder.create(maxTaxaNumber,barcode.getTaxaIndex()));
                }
                else if(taxaDistribution.totalDepth()==1) { //need to go from singleton to expandable TaxaDistribution
                    masterTagTaxaMap.put(tag,TaxaDistBuilder.create(taxaDistribution).increment(barcode.getTaxaIndex()));
                } else {
                    taxaDistribution.increment(barcode.getTaxaIndex());
                }
                if (allReads % 1000000 == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read");
                }
            }
            myLogger.info("Total number of reads in lane=" + allReads);
            myLogger.info("Total number of good barcoded reads=" + goodBarcodedReads);
            myLogger.info("Timing process (sorting, collapsing, and writing TagCount to file).");
            System.out.println("tagCntMap size"+masterTagTaxaMap.size());
            myLogger.info("Process took " + (System.nanoTime() - time)/1e6 + " milliseconds.");
            br.close();
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);
            e.printStackTrace();
        }
    }

    /**
     * Method for reading FastQ four line structure, and returning a string array with [sequence, qualityScore]
     */
    static String[] readFastQBlock(BufferedReader bw, int currentRead) throws IOException {
        //consider converting this into a stream of String[]
        String[] result=new String[2];
        try{
            bw.readLine();
            result[0]=bw.readLine();
            bw.readLine();
            result[1]=bw.readLine();
            if(result[0]==null) {
                return null;
            }
            return result;
        } catch (IOException e) {
            e.printStackTrace();
            myLogger.error("Unable to correctly parse the sequence and quality score near line: " + currentRead*4
                    + " from fastq file.  Your fastq file may have been corrupted.");
            return null;
        }
    }

    private void removeSecondCutSitesFromMap(GBSEnzyme enzyme) {
        //this is a little tricky as you cannot add entries at the same time as removing entries to a map
        System.out.println("DiscoveryTBTPlugin.removeSecondCutSitesFromMap started Initial Size:"+tagCntMap.size());
        String[] likelyReadEnd=enzyme.likelyReadEnd();
        Map<Tag,TaxaDistribution> shortTags=new HashMap<>(tagCntMap.size()/5);
        int belowMinSize=0, shortExisting=0;
        for (Tag origTag : tagCntMap.keySet()) {
            int minCutSite=Integer.MAX_VALUE;
            for (String potentialCutSite : likelyReadEnd) {
                int p = origTag.sequence().indexOf(potentialCutSite, 1);
                if(p>0) minCutSite=Math.min(minCutSite,p);
            }
            if (minCutSite!=Integer.MAX_VALUE && minCutSite > 1) {
                if(minCutSite<minimumTagLength()) {
                    tagCntMap.remove(origTag);
                    belowMinSize++;
                    continue;
                }
                TaxaDistribution currentTaxaDist=tagCntMap.remove(origTag);
                Tag t = TagBuilder.instance(origTag.sequence().substring(0, minCutSite + enzyme.readEndCutSiteRemnantLength())).build();
                TaxaDistribution existingTD=shortTags.get(t);
                if(existingTD!=null) {
                    if(currentTaxaDist==null || existingTD==null) {
                        System.out.println("We have a problem");
                    }
                    currentTaxaDist=TaxaDistBuilder.combine(currentTaxaDist,existingTD);
                    shortExisting++; //System.out.println("existingTD="+t);
                }
                shortTags.put(t, currentTaxaDist);
                //if(ttemp!=null) shortExisting++;
               // System.out.println(origTag+" -> "+t);
            }
        }
        System.out.println("After removal tagCntMap.size() = " + tagCntMap.size());
        System.out.println("After removal shortTags.size() = " + shortTags.size());
        System.out.println("belowMinSize = " + belowMinSize);
        System.out.println("shortExisting = " + shortExisting);
        tagCntMap.putAll(shortTags);
        System.out.println("After combining again tagCntMap.size() = " + tagCntMap.size());
    }

// The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(DiscoveryTBTPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
//    public <Type> runPlugin(DataSet input) {
//        return (<Type>) performFunction(input).getData(0).getData();
//    }

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
     * Maximum Tag Length
     *
     * @return Maximum Tag Length
     */
    public Integer maximumTagLength() {
        return myMaxTagLength.value();
    }

    /**
     * Set Maximum Tag Length. Maximum Tag Length
     *
     * @param value Maximum Tag Length
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin maximumTagLength(Integer value) {
        myMaxTagLength = new PluginParameter<>(myMaxTagLength, value);
        return this;
    }

    /**
     * Minimum Tag Length
     *
     * @return Minimum Tag Length
     */
    public Integer minimumTagLength() {
        return myMinTagLength.value();
    }

    /**
     * Set Minimum Tag Length. Minimum Tag Length
     *
     * @param value Minimum Tag Length
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin minimumTagLength(Integer value) {
        myMinTagLength = new PluginParameter<>(myMinTagLength, value);
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
     * Output Database File
     *
     * @return Output Database File
     */
    public String outputDatabaseFile() {
        return myOutputDB.value();
    }

    /**
     * Set Output Database File. Output Database File
     *
     * @param value Output Database File
     *
     * @return this plugin
     */
    public DiscoveryTBTPlugin outputDatabaseFile(String value) {
        myOutputDB = new PluginParameter<>(myOutputDB, value);
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
     *
     * @return Maximum Map Memory in Mb
     */
    public Integer maximumMapMemoryInMb() {
        return myMaxMapMemoryInMb.value();
    }

    /**
     * Set Maximum Map Memory in Mb. Maximum size for the
     * tag distribution map in Mb
     *
     * @param value Maximum Map Memory in Mb
     *
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

    /**
     * This ConcurrentHashMap constrain the size of the map, and purges low distribution count tags when the size needs
     * to be reduced.
     */
    static class TagDistributionMap extends ConcurrentHashMap<Tag,TaxaDistribution> {
        private final long maxMemorySize;
        private int minDepthToRetainInMap=1;
        private LongAdder putCntSinceMemoryCheck=new LongAdder();  //since we don't need an exact count of the puts,
        //this may be overkill
        private long checkFreq;  //numbers of puts to check size

        TagDistributionMap(int initialCapacity, long maxMemorySize) {
            super(initialCapacity);
            this.maxMemorySize=maxMemorySize;
            checkFreq=(maxMemorySize/(100L*10L));  //this is checking roughly to keep the map within 10% of the max
        }

        @Override
        public TaxaDistribution put(Tag key, TaxaDistribution value) {
            putCntSinceMemoryCheck.increment();
            if(putCntSinceMemoryCheck.longValue()>checkFreq) {
                checkMemoryAndReduceIfNeeded();
                putCntSinceMemoryCheck.reset();
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
            entrySet().parallelStream()
                    .filter(e -> e.getValue().totalDepth()<minCnt)
                    .forEach(e -> remove(e.getKey()));
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
}


