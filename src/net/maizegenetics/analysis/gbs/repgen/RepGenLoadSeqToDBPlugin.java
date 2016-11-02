package net.maizegenetics.analysis.gbs.repgen;

import com.google.common.collect.Multimap;
import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
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
import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Develops a discovery rGBS database based on a folder of sequencing files
 *
 * Keeps only good reads having no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 * 
 * Originally the reference throughout was to "tag". This is being changed
 * to "kmer" as the pipeline is a kmer alignment process.  
 *
 * @author Ed Buckler
 */
public class RepGenLoadSeqToDBPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(RepGenLoadSeqToDBPlugin.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing FASTQ files in text or gzipped text.\n"
                    + "     NOTE: Directory will be searched recursively and should\n"
                    + "     be written WITHOUT a slash after its name.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<Integer> myKmerLength = new PluginParameter.Builder<>("kmerLength", 150, Integer.class).guiName("Maximum Kmer Length")
            .description("Specified length for each kmer to process").build();
    private PluginParameter<Integer> myMinKmerLength = new PluginParameter.Builder<>("minKmerL", 120, Integer.class).guiName("Minimum Kmer Length")
            .description("Minimum kmer Length after second cut site is removed").build();
    private PluginParameter<Integer> myMinKmerCount = new PluginParameter.Builder<>("c", 10, Integer.class).guiName("Min Kmer Count")
            .description("Minimum kmer count").build();
    private PluginParameter<String> myOutputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Output Database File").required(true).outFile()
            .description("Output Database File").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();
    private PluginParameter<Integer> myMaxKmerNumber = new PluginParameter.Builder<>("mxKmerNum", 50000000, Integer.class).guiName("Maximum Kmer Number").required(false)
            .description("Maximum number of kmers").build();
    private PluginParameter<Integer> myBatchSize = new PluginParameter.Builder<>("batchSize", 8, Integer.class).guiName("Batch size of fastq files").required(false)
            .description("Number of flow cells being processed simultaneously").build();
    private PluginParameter<Boolean> myDeleteOldData = new PluginParameter.Builder<Boolean>("deleteOldData",false,Boolean.class).guiName("Delete Old Data")
            .description("Delete existing SNP quality data from db tables").build();
    LongAdder roughTagCnt = new LongAdder();

    private TagDistributionMap tagCntMap;
    private boolean taglenException;
    protected static int readEndCutSiteRemnantLength;
    String[] likelyReadEndStrings;

    public RepGenLoadSeqToDBPlugin() {
        super(null, false);
    }

    public RepGenLoadSeqToDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    private long[] calcTagMapStats(TagDistributionMap tagCntMap) {
        long totalDepth=0, memory=0;
        int cnt=0;
        for (Map.Entry<Tag, TaxaDistribution> entry : tagCntMap.entrySet()) {
            memory+=entry.getValue().memorySize();
            memory+=25;
            totalDepth+=entry.getValue().totalDepth();
            cnt++;
        }
        int currentSize = tagCntMap.size();
        memory+=tagCntMap.size()*2*16;  //estimate for the map size
        long[] stats={currentSize,memory, totalDepth,totalDepth/cnt};
        System.out.printf("Map Tags:%,d  Memory:%,d  TotalDepth:%,d  AvgDepthPerTag:%d%n",stats[0],stats[1],stats[2],stats[3]);
        return stats;
    }
    @Override
    public void postProcessParameters() {

    }
    @Override
    public DataSet processData(DataSet input) {
        int batchSize = myBatchSize.value();
        float loadFactor = 0.95f;
        tagCntMap = new TagDistributionMap (myMaxKmerNumber.value(),loadFactor, 128, this.minKmerCount());
        try {
            TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSUtils.sampleNameField, new HashMap<>(), true);
            Map<String,Taxon> fileTaxaMap=TaxaListIOUtils.getUniqueMapOfTaxonByAnnotation(masterTaxaList,GBSUtils.fileNameField)
                    .orElseThrow(() -> new IllegalArgumentException("Error: Same file points more than one taxon in the KeyFile"));
            //Get the list of fastq files
            List<Path> directoryFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
            if(directoryFiles.isEmpty()) {
                myLogger.warn("No files matching:"+ GBSUtils.inputFileGlob);
                return null;
            }
            // Cull files that are not represented in the given key file
            List<Path> inputSeqFiles =  directoryFiles.stream()
                    .peek(path -> System.out.println(path.getFileName().toString()))
                    .filter(path -> fileTaxaMap.containsKey(path.getFileName().toString()))
                    .collect(Collectors.toList());
            if (inputSeqFiles.size() == 0) {
                System.out.println("RepGenLoadSeqToDB:processData - found NO files represented in key file.");
                System.out.println("Please verify your file names are formatted correctly and that your key file contains the required headers.");
                return null; // no files in this directory to process
            }
            System.out.println("Found " + inputSeqFiles.size() + " files to process");
            //Files in a batch have roughly the same size
            //inputSeqFiles = this.sortFastqBySize(inputSeqFiles);
            int batchNum = inputSeqFiles.size()/batchSize;
            if (inputSeqFiles.size()%batchSize != 0) batchNum++;

            // Check if user wants to clear existing db.
            RepGenDataWriter tdw = null;
            if (Files.exists(Paths.get(myOutputDB.value()))) {
                if (deleteOldData()) {
                    try {
                        Files.delete(Paths.get(myOutputDB.value()));
                    } catch (Exception exc){
                        System.out.println("Error when trying to delete database file: " + myOutputDB.value());
                        System.out.println("File delete error: " + exc.getMessage());
                        return null;
                    }
                } else {
                 // We'll append data to new DB if all new taxa are contained in db
                    tdw=new RepGenSQLite(myOutputDB.value());
                    TaxaList oldTaxaList = tdw.getTaxaList();
                    boolean sameMasterTaxaList = true;
                    Taxon badTaxon = null;
                    for (Taxon taxa : masterTaxaList) {
                        if (!oldTaxaList.contains(taxa)) {
                            sameMasterTaxaList = false;
                            badTaxon = taxa;
                            break;
                        }
                    }
                    if (!sameMasterTaxaList) {
                        myLogger.error("Taxa list in keyfile contains taxa not in db:  " + badTaxon +
                            ".\nEither delete existing DB, use a new DB output file, or use keyfile with entries matching existing taxa list.\n");
                        ((TagDataSQLite)tdw).close();
                        return null;
                    }
                    // Grab existing data from db, append to empty tagCntMap
                    Map<Tag, TaxaDistribution> existingTDM = tdw.getAllTagsTaxaMap();
                    tagCntMap.putAll(existingTDM);
                    tdw.clearTagTaxaDistributionData(); // clear old data - it will be re-added at the end.
                }
            }
            if (tdw == null) tdw=new RepGenSQLite(myOutputDB.value());
            taglenException = false;
            for (int i = 0; i < inputSeqFiles.size(); i+=batchSize) {
                int end = i+batchSize;
                if (end > inputSeqFiles.size()) end = inputSeqFiles.size();
                ArrayList<Path> sub = new ArrayList();
                for (int j = i; j < end; j++) sub.add(inputSeqFiles.get(j));
                System.out.println("\nStart processing batch " + String.valueOf(i/batchSize+1));
                sub//.parallelStream()
                    .forEach(inputSeqFile -> {
                        try {
                            //processFastQFile(masterTaxaList,keyPath, inputSeqFile, minimumQualityScore(), tagCntMap, kmerLength());
                            int taxaIndex=masterTaxaList.indexOf(fileTaxaMap.get(inputSeqFile.getFileName().toString()));
                            processFastQ(inputSeqFile, taxaIndex, masterTaxaList, tagCntMap, kmerLength(), minimumQualityScore());
                        } catch (StringIndexOutOfBoundsException oobe) {
                            oobe.printStackTrace();
                            myLogger.error(oobe.getMessage());
                            setTagLenException();
                            return;
                        }
                });
                if (taglenException == true) return null; // Tag length failure from processFastQ - halt processing

                System.out.println("\nKmers are added from batch "+String.valueOf(i/batchSize+1) + ". Total batch number: " + batchNum);
                int currentSize = tagCntMap.size();
                System.out.println("Current number: " + String.valueOf(currentSize) + ". Max kmer number: " + String.valueOf(myMaxKmerNumber.value()));
                System.out.println(String.valueOf((float)currentSize/(float)myMaxKmerNumber.value()) + " of max tag number");

                if (currentSize > 0) { // calcTagMapStats() gets "divide by 0" error when size == 0
                    this.calcTagMapStats(tagCntMap);
                    System.out.println();
                    //make sure don't lose rare ones, need to set maxTagNumber large enough
                    removeTagsWithoutReplication(tagCntMap);
                    if (tagCntMap.size() == 0) {
                        System.out.println("WARNING:  After removing tags without replication, there are NO  tags left in the database");
                    } else {
                        this.calcTagMapStats(tagCntMap);
                        System.out.println();
                        System.out.println("Kmer number is reduced to " + tagCntMap.size()+"\n");
                    }
                    this.roughTagCnt.reset();
                    this.roughTagCnt.add(tagCntMap.size());
                } else {
                    System.out.println("WARNING: Current tagcntmap size is 0 after processing batch " + String.valueOf(i/batchSize+1) );
                }

                System.out.println("Total memory: "+ String.valueOf((double)(Runtime.getRuntime().totalMemory()/1024/1024/1024))+" Gb");
                System.out.println("Free memory: "+ String.valueOf((double)(Runtime.getRuntime().freeMemory()/1024/1024/1024))+" Gb");
                System.out.println("Max memory: "+ String.valueOf((double)(Runtime.getRuntime().maxMemory()/1024/1024/1024))+" Gb");
                System.out.println("\n");
            }
            System.out.println("\nAll the batch are processed");
            tagCntMap.removeTagByCount(myMinKmerCount.value());
            System.out.println("By removing kmers with minCount of " + myMinKmerCount.value() + "Kmer number is reduced to " + tagCntMap.size()+"\n");

            tdw.putTaxaList(masterTaxaList);
            tdw.putAllTag(tagCntMap.keySet());
            tdw.putTaxaDistribution(tagCntMap);
            ((RepGenSQLite)tdw).close();  //todo autocloseable should do this but it is not working.
        } catch(Exception e) {
            e.printStackTrace();
        }
        return new DataSet(new Datum("TagMap",tagCntMap,""),this);
    }

    private void processFastQ(Path fastqFile, int taxaIndex, TaxaList masterTaxaList,
                              TagDistributionMap masterTagTaxaMap, int preferredTagLength, int minQual) throws StringIndexOutOfBoundsException{
        int allReads=0, goodBarcodedReads = 0, lowQualityReads = 0;
        int maxTaxaNumber=masterTaxaList.size();
        int checkSize = 10000000;
        myLogger.info("processing file " + fastqFile.toString());
        try {
            int qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;

            while ((seqAndQual=GBSUtils.readFastQBlock(br,allReads)) != null) {
                allReads++;
                // Check for and toss low quality reads
                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual, qualityScoreBase)<(preferredTagLength)){
                    	lowQualityReads++;
                    	continue;
                    }
                }

                if (seqAndQual[0].length() < preferredTagLength) {
                	String errMsg = "\n\nERROR processing " + fastqFile.toString() + "\n" +
                			"Reading entry number " + allReads + " fails the length test.\n" +
                			"Sequence length " + seqAndQual[0].length() +
                			" is less then maxKmerLength " + preferredTagLength + ".\n" +
                			"Re-run your files with either a shorter mxKmerL value or a higher minimum quality score.\n";
                	throw new StringIndexOutOfBoundsException(errMsg);
                }
                
                Tag tag = TagBuilder.instance(seqAndQual[0].substring(0,preferredTagLength)).build();

                if(tag==null) continue;   //null occurs when any base was not A, C, G, T
                goodBarcodedReads++;
                TaxaDistribution taxaDistribution=masterTagTaxaMap.get(tag);
                if(taxaDistribution==null) {
                    masterTagTaxaMap.put(tag,TaxaDistBuilder.create(maxTaxaNumber,taxaIndex));
                    this.roughTagCnt.increment();
                } else {
                    taxaDistribution.increment(taxaIndex);
                }
                if (allReads % checkSize == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read. Current tag count:" + this.roughTagCnt);
                }
            }
            myLogger.info("Summary for "+fastqFile.toString()+"\n"+
                    "Total number of reads in lane=" + allReads +"\n"+
                    "Total number of good barcoded reads=" + goodBarcodedReads+"\n"+
                    "Total number of low quality reads=" + lowQualityReads+"\n"+
                    "Timing process (sorting, collapsing, and writing TagCount to file)."+"\n"+
                    "Process took " + (System.nanoTime() - time)/1e6 + " milliseconds.");
            System.out.println("tagCntMap size: "+masterTagTaxaMap.size());
            br.close();
        } catch (StringIndexOutOfBoundsException oobe) {
        	throw oobe; // pass it up to print error and stop processing
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);

            e.printStackTrace();
        }
    }


    /**
     * This method removes all tags are are never repeated in a single sample (taxa).  The concept is that
     * all biologically real tag should show up twice somewhere.  This could be called at the end of every
     * flowcell to test all the novel tags.
     */
    private static void removeTagsWithoutReplication (TagDistributionMap masterTagTaxaMap) {
        int currentSize = masterTagTaxaMap.size();
        int minTaxa=2;
        System.out.println("Starting removeTagsWithoutReplication. Current tag number: " + currentSize);
        LongAdder tagsRemoved=new LongAdder();
        masterTagTaxaMap.entrySet().parallelStream().forEach(t -> {
            TaxaDistribution td = t.getValue();
            if(td.totalDepth()<2*minTaxa) {
                masterTagTaxaMap.remove(t.getKey());
                tagsRemoved.increment();
            }
            else if(IntStream.of(td.depths()).filter(depth -> depth>1).count()<minTaxa) {
                masterTagTaxaMap.remove(t.getKey());
                tagsRemoved.increment();
            }
        });
        System.out.println("Finished removeTagsWithoutReplication.  tagsRemoved = " + tagsRemoved + ". Current tag number: " + String.valueOf(currentSize-tagsRemoved.intValue()));
    }

    public void setTagLenException() {
    	taglenException = true;
    }

// The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GBSSeqToTagDBPlugin.class);
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
    public RepGenLoadSeqToDBPlugin inputDirectory(String value) {
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
    public RepGenLoadSeqToDBPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }


    /**
     * Maximum Tag Length
     *
     * @return Maximum Tag Length
     */
    public Integer kmerLength() {
        return myKmerLength.value();
    }

    /**
     * Set Maximum Tag Length. Maximum Tag Length
     *
     * @param value Maximum Tag Length
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin kmerLength(Integer value) {
        myKmerLength = new PluginParameter<>(myKmerLength, value);
        return this;
    }

    /**
     * Minimum Tag Length
     *
     * @return Minimum Tag Length
     */
    public Integer minimumKmerLength() {
        return myMinKmerLength.value();
    }

    /**
     * Set Minimum Tag Length. Minimum Tag Length
     *
     * @param value Minimum Tag Length
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin minimumKmerLength(Integer value) {
        myMinKmerLength = new PluginParameter<>(myMinKmerLength, value);
        return this;
    }

    /**
     * Minimum tag count
     *
     * @return Min Tag Count
     */
    public Integer minKmerCount() {
        return myMinKmerCount.value();
    }

    /**
     * Set Min Tag Count. Minimum tag count
     *
     * @param value Min Tag Count
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin minKmerCount(Integer value) {
        myMinKmerCount = new PluginParameter<>(myMinKmerCount, value);
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
    public RepGenLoadSeqToDBPlugin outputDatabaseFile(String value) {
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
    public RepGenLoadSeqToDBPlugin minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }

    /**
     * Set maximum number of kmers
     * @param value
     * @return
     */
    public RepGenLoadSeqToDBPlugin maximumKmerNumber(Integer value) {
        myMaxKmerNumber = new PluginParameter<>(myMaxKmerNumber, value);
        return this;
    }

    /**
     * Set number of Fastq files processed simultaneously
     * @param value
     * @return
     */
    public RepGenLoadSeqToDBPlugin batchSize(Integer value) {
        myBatchSize = new PluginParameter<>(myBatchSize, value);
        return this;
    }
    /**
     * Delete exisiting  DB
     *
     * @return deleteOldData
     */
    public Boolean deleteOldData() {
        return myDeleteOldData.value();
    }

    /**
     * Set Delete old data flag.  True indicates we want the
     * db tables cleared
     *
     * @param value true/false - whether to delete data
     *
     * @return this plugin
     */
    public RepGenLoadSeqToDBPlugin deleteOldData(Boolean value) {
        myDeleteOldData = new PluginParameter<>(myDeleteOldData, value);
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
        private final int maxTagNum;
        private int minDepthToRetainInMap=2;
        private final int minCount;

        TagDistributionMap (int maxTagNumber, float loadFactor, int concurrencyLevel, int minCount) {
            super((maxTagNumber*2), loadFactor, concurrencyLevel);
            maxTagNum = maxTagNumber;
            this.minCount = minCount;
        }

        @Override
        public TaxaDistribution put(Tag key, TaxaDistribution value) {
            return super.put(key, value);
        }

        public synchronized void removeTagByCount(int minCnt) {
            entrySet().parallelStream()
                    .filter(e -> e.getValue().totalDepth()<minCnt)
                    .forEach(e -> remove(e.getKey()));
        }


        public long estimateMapMemorySize() {
            long size=0;
            int cnt=0;
            for (Entry<Tag, TaxaDistribution> entry : entrySet()) {
                size+=8+16+1; //Tag size
                size+=16; //Map references
                size+=entry.getValue().memorySize();
                cnt++;
                if(cnt>10000) break;
        }
            long estSize=(size()/cnt)*size;
            return estSize;
        }

        public long[] depthDistribution() {
            long[] base2bins=new long[34];
            int cnt=0;
            for (Entry<Tag, TaxaDistribution> entry : entrySet()) {
                base2bins[31-Integer.numberOfLeadingZeros(entry.getValue().totalDepth())]++;
                cnt++;
               // if(cnt>100000) break;
            }
            return base2bins;
        }
    }
}


