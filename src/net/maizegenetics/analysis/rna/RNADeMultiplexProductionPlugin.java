package net.maizegenetics.analysis.rna;

import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import net.maizegenetics.analysis.gbs.v2.GBSEnzyme;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.analysis.gbs.v2.ProductionSNPCallerPluginV2;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.TaxaTissueDist;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;


public class RNADeMultiplexProductionPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(ProductionSNPCallerPluginV2.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing fastq AND/OR qseq files.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myInputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Input GBS Database").required(true).inFile()
            .description("Input Database file if using SQLite").build();
    private PluginParameter<FindMatchByKmers.MatchType> myMatchingType = new PluginParameter.Builder<>("matchType", null, FindMatchByKmers.MatchType.class)
            .guiName("Matching approach").required(true)
            .range(Range.encloseAll(Arrays.asList(FindMatchByKmers.MatchType.values())))
            .description("Approach used for matching the reads to the contigs").build();
    private PluginParameter<Integer> myKmerForMatching = new PluginParameter.Builder<>("kmer", 16, Integer.class).guiName("Minimum Kmer Length")
            .description("Minimum length for kmer to be kept").build();
    private PluginParameter<Integer> myBatchSize = new PluginParameter.Builder<>("batchSize", 8, Integer.class).guiName("Batch size of fastq files").required(false)
            .description("Number of flow cells being processed simultaneously").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();


    private TagDataWriter tdw = null;
    private Map<String,int[]> taxatissueCntMap;
    private TagTissueDistributionMap tagTaxaTissueMap;
    private FindMatchByKmers findMatchByKmers;
 
    protected static int readEndCutSiteRemnantLength;
    private static String myEnzyme="ignore";
    private static int maxTaxa;
    private static int maxTissue;

    //Documentation of read depth per sample (one recorded per replicate)
    // Treemap is synchronized as multiple threads may increment values.
    private Map<String, Integer> rawReadCountsMap = new TreeMap<>();
    private Map<String, Integer> rawReadCountsForFullSampleName = Collections.synchronizedMap(rawReadCountsMap);
    private Map<String, Integer> matchedReadCountsMap = new TreeMap<>();
    private Map<String, Integer> matchedReadCountsForFullSampleName = Collections.synchronizedMap(matchedReadCountsMap);

    private boolean taglenException;
    
    public RNADeMultiplexProductionPlugin() {
        super(null, false);
    }

    public RNADeMultiplexProductionPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        GBSEnzyme enzyme = new GBSEnzyme(myEnzyme);
        readEndCutSiteRemnantLength = enzyme.readEndCutSiteRemnantLength();
    }

    @Override
    public DataSet processData(DataSet input) {
        int batchSize = batchSize();
        Path keyPath= Paths.get(keyFile()).toAbsolutePath();
        List<Path> directoryFiles= DirectoryCrawler.listPaths(GBSUtils.inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
        if(directoryFiles.isEmpty()) {
            myLogger.warn("No files matching:"+GBSUtils.inputFileGlob);
            return null;
        }
        List<Path> inputSeqFiles =directoryFiles;
//        List<Path> inputSeqFiles = GBSUtils.culledFiles(directoryFiles,keyPath);
//        if (inputSeqFiles.size() == 0) return null; // no files to process

        tdw =new TagDataSQLite(myInputDB.value());



        TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSUtils.sampleNameField, new HashMap<>(), true);
        Multimap<String,Taxon> temp=TaxaListIOUtils.getMapOfTaxonByAnnotation(masterTaxaList, "FileName");
        Map<String,Taxon> fileToTaxonMap=new HashMap<>();
        temp.entries().stream().forEach(entries -> {
            Taxon t=fileToTaxonMap.put(entries.getKey(), entries.getValue());
            if(t!=null) System.err.println("Duplicate file");
        });  //wrote over duplicated

        ArrayList<String> masterTissueList = (ArrayList)TaxaListIOUtils.readTissueAnnotationFile(keyFile(), GBSUtils.tissueNameField);
        // Write tissue list to the db
        tdw.putTaxaList(masterTaxaList);
        tdw.putAllTissue(masterTissueList);

        maxTaxa = masterTaxaList.size();
        maxTissue = masterTissueList.size();

        //get all the possible tags (=sequences = contigs)
        findMatchByKmers=new FindMatchByKmers(tdw.getTags(), matchingType(), kmerForMatching(), 10);
        taxatissueCntMap=new ConcurrentHashMap<>(maxTissue*maxTaxa*2);
               
        writeInitialTaxaReadCounts(masterTaxaList); // initialize synchronized maps

        int batchNum = inputSeqFiles.size()/batchSize;
       
        if (inputSeqFiles.size() % batchSize !=0) batchNum++;
        System.out.println("ProductionSNPCallerPluginV2: Total batches to process: " + batchNum);


        taglenException = false;
        for (int idx = 0; idx < inputSeqFiles.size(); idx+=batchSize) {
            int end = idx+batchSize;
            if (end > inputSeqFiles.size()) end = inputSeqFiles.size();
            ArrayList<Path> sub = new ArrayList<Path>();
            for (int jdx = idx; jdx < end; jdx++) sub.add(inputSeqFiles.get(jdx));
            System.out.println("\nStart processing batch " + String.valueOf(idx/batchSize+1));
            sub.parallelStream()
            .forEach(inputSeqFile -> {
                try {
                    processFastQ(inputSeqFile,fileToTaxonMap.get(inputSeqFile.getFileName().toString()),minQualScore());

                } catch (StringIndexOutOfBoundsException oobe) {
                    oobe.printStackTrace();
                    myLogger.error(oobe.getMessage());
                    setTagLenException();
                    return;
                }              
            });
            if (taglenException == true) return null; // Tag length failure from processFastQ - halt processing
         
            System.out.println("\nFinished processing batch " + String.valueOf(idx/batchSize+1));
        }
        
        // put counts to database.
        taxatissueCntMap.entrySet().stream().forEach(entry -> {
            String[] taxaTissueNames = entry.getKey().split(",");
            tdw.putTaxaTissueDistribution(taxaTissueNames[0],taxaTissueNames[1],findMatchByKmers.getTags(),entry.getValue());
        });

//        writeReadsPerSampleReports(tagTaxaTissueMap.size());
        return new DataSet(new Datum("DB Name",inputDB(),""),null);
    }


    private void processFastQ(Path fastqFile, Taxon taxon, int minQual) throws StringIndexOutOfBoundsException {
        String tissuetaxaKey=taxon.getName()+
                ","+taxon.getAnnotation().getTextAnnotation("Tissue")[0];
        taxatissueCntMap.putIfAbsent(tissuetaxaKey, new int[findMatchByKmers.totalNumberOfTags()]);
        int allReads=0, goodBarcodedReads = 0, lowQualityReads = 0;
        int shortReads = 0;
        try {
            int qualityScoreBase=GBSUtils.determineQualityScoreBase(fastqFile);
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;
            String likelyReadEnd = "AGATCGGA";
            while ((seqAndQual=GBSUtils.readDeMultiPlexFastQBlock(br, allReads)) != null) {
                allReads++;
                int adapterStart = 0;
                // Decode barcode using the current sequence & quality  score

                // Find barcode for smallRNA (taxa comes from barcode)
                // Strip off the barcode - original sequence did NOT have it,
                // this was added when we read the fastQ block.
                // We to revert the sequence to what is in the fastQ file
                // so the quality scores match up correctly.
                String sequence = seqAndQual[0];
                // Substring it once more to remove the common adaptor
                // and everything beyond it

                adapterStart = sequence.indexOf(likelyReadEnd);
                if (adapterStart > 0) {
                    sequence = sequence.substring(0, adapterStart-1);
                }
//                if (sequence.length() < minKmerLen){
//                    System.out.println("LCJ - found short read, seq: " + sequence);
//                     shortReads++;
//                    continue;
//                }
                if(minQual>0) {
                    //todo move getFirstLowQualityPos into this class?
                    if(BaseEncoder.getFirstLowQualityPos(seqAndQual[1],minQual, qualityScoreBase) < sequence.length()){
                        lowQualityReads++;
                        continue;
                    }
                }
                OptionalInt hitIndex=findMatchByKmers.getMatchIndex(sequence);
                if(hitIndex.isPresent()) {
                    goodBarcodedReads++;
                    taxatissueCntMap.get(tissuetaxaKey)[hitIndex.getAsInt()]++;
                }

                if (allReads % 1000000 == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read");
                }
            }
            myLogger.info("Total number of reads in lane=" + allReads);
            myLogger.info("Total number of good barcoded reads=" + goodBarcodedReads);
            myLogger.info("Total number of low quality reads=" + lowQualityReads);
            myLogger.info("Timing process (sorting, collapsing, and writing TagCount to file).");
            myLogger.info("Process took " + (System.nanoTime() - time)/1e6 + " milliseconds for file " + fastqFile.toString());
            br.close();
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);
            e.printStackTrace();
        }
    }



//    private void writeReadsPerSampleReports(int tagsProcessed) {
//        myLogger.info("\nWriting ReadsPerSample log file...");
//        String outFileS = myOutputDir + File.separator + (new File(keyFile())).getName();
//        outFileS = outFileS.replaceAll(".txt", "_ReadsPerSample.log");
//        outFileS = outFileS.replaceAll("_key", "");
//        try {
//                String msg = "ReadsPerSample log file: " + outFileS;
//                myLogger.info(msg);
//            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFileS))), 65536);
//            bw.write("FullSampleName\t\t\tgoodBarcodedReads\tgoodReadsMatchedToDataBase\n");
//            for (String fullSampleName : rawReadCountsForFullSampleName.keySet()) {
//                bw.write(fullSampleName + "\t" + rawReadCountsForFullSampleName.get(fullSampleName) + "\t\t" + matchedReadCountsForFullSampleName.get(fullSampleName) + "\n");
//            }
//            bw.close();
//        } catch (Exception e) {
//            myLogger.error("Couldn't write to ReadsPerSample log file: " + e);
//            e.printStackTrace();
//            System.exit(1);
//        }
//        myLogger.info("\n\nTotal number of SNPs processed with minimum quality score " + minimumQualityScore() + " was " + tagsProcessed + ".\n");
//        myLogger.info("   ...done\n");
//    }
    
    private void writeInitialTaxaReadCounts(TaxaList tl) {
        tl.stream() // Add initial taxa names with count of 0 to synchronized maps
        .forEach(taxon -> {
                 rawReadCountsForFullSampleName.put(taxon.getName(), 0); 
             matchedReadCountsForFullSampleName.put(taxon.getName(), 0);
        });
    }

    public void setTagLenException() {
        taglenException = true;
    }
    

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Production SNP Caller";
    }

    @Override
    public String getToolTipText() {
        return "Production SNP Caller";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(RNADeMultiplexProductionPlugin.class);
     }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
    }

    /**
     * Input directory containing fastq AND/OR qseq files.
     *
     * @return Input Directory
     */
    public String inputDir() {
        return myInputDir.value();
    }

    /**
     * Set Input Directory. Input directory containing fastq
     * AND/OR qseq files.
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin inputDir(String value) {
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
    public RNADeMultiplexProductionPlugin keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }

    /**
     * Input Database file if using SQLite
     *
     * @return Input GBS Database
     */
    public String inputDB() {
        return myInputDB.value();
    }

    /**
     * Set Input GBS Database. Input Database file if using
     * SQLite
     *
     * @param value Input GBS Database
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin inputDB(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
        return this;
    }

    /**
     * Approach used for matching the reads to the contigs
     *
     * @return Matching approach
     */
    public FindMatchByKmers.MatchType matchingType() {
        return myMatchingType.value();
    }

    /**
     * Set Matching approach. Approach used for matching the
     * reads to the contigs
     *
     * @param value Matching approach
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin matchingType(FindMatchByKmers.MatchType value) {
        myMatchingType = new PluginParameter<>(myMatchingType, value);
        return this;
    }

    /**
     * Minimum length for kmer to be kept
     *
     * @return Minimum Kmer Length
     */
    public Integer kmerForMatching() {
        return myKmerForMatching.value();
    }

    /**
     * Set Minimum Kmer Length. Minimum length for kmer to
     * be kept
     *
     * @param value Minimum Kmer Length
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin kmerForMatching(Integer value) {
        myKmerForMatching = new PluginParameter<>(myKmerForMatching, value);
        return this;
    }

    /**
     * Number of flow cells being processed simultaneously
     *
     * @return Batch size of fastq files
     */
    public Integer batchSize() {
        return myBatchSize.value();
    }

    /**
     * Set Batch size of fastq files. Number of flow cells
     * being processed simultaneously
     *
     * @param value Batch size of fastq files
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin batchSize(Integer value) {
        myBatchSize = new PluginParameter<>(myBatchSize, value);
        return this;
    }

    /**
     * Minimum quality score within the barcode and read length
     * to be accepted
     *
     * @return Minimum quality score
     */
    public Integer minQualScore() {
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
    public RNADeMultiplexProductionPlugin minQualScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }

    /**
     * Set Minimum quality score. Minimum quality score within
     * the barcode and read length to be accepted
     *
     * @param value Minimum quality score
     *
     * @return this plugin
     */
    public RNADeMultiplexProductionPlugin minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }
    
    static class TagTissueDistributionMap extends ConcurrentHashMap<Tag,TaxaTissueDist> {
        private final int maxTagNum;
        private int minDepthToRetainInMap=2;
        private final int minCount; // minimum count before tag is tossed
        
        TagTissueDistributionMap (int maxTagNumber, float loadFactor, int concurrencyLevel, int minCount) {
            super((maxTagNumber*2), loadFactor, concurrencyLevel);
            maxTagNum = maxTagNumber;
            this.minCount = minCount;
        }
        
        // This is putting the value into the map, not writing the
        // counts.  But we need to be incrementing concurrently -
        // where do we "put", where do we "increment".
        @Override
        public TaxaTissueDist put(Tag key, TaxaTissueDist value) {
            return super.put(key, value);
        }
    }

}
