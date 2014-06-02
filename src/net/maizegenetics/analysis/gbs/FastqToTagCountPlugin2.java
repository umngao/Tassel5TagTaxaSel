package net.maizegenetics.analysis.gbs;

import com.google.common.collect.BiMap;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TaxaDist;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Derives a tagCount list for each fastq file in the input directory.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 *
 */
public class FastqToTagCountPlugin2 extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FastqToTagCountPlugin2.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<String>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing FASTQ files in text or gzipped text.\n"
                    + "     NOTE: Directory will be searched recursively and should\n"
                    + "     be written WITHOUT a slash after its name.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<String>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myEnzyme = new PluginParameter.Builder<String>("e", null, String.class).guiName("Enzyme").required(true)
            .description("Enzyme used to create the GBS library, if it differs from the one listed in the key file").build();
    private PluginParameter<Integer> myMaxGoodReads = new PluginParameter.Builder<Integer>("s", 300000000, Integer.class).guiName("Max Good Reads")
            .description("Max good reads per lane").build();
    private PluginParameter<Integer> myMinTagCount = new PluginParameter.Builder<Integer>("c", 1, Integer.class).guiName("Min Tag Count")
            .description("Minimum tag count").build();
    private PluginParameter<String> myOutputDir = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output Directory").required(true).outDir()
            .description("Output directory to contain .cnt files (one per FASTQ file (one per FASTQ file)").build();
    private PluginParameter<Integer> myMinQualScore = new PluginParameter.Builder<>("mnQS", 0, Integer.class).guiName("Minimum quality score").required(false)
            .description("Minimum quality score within the barcode and read length to be accepted").build();

    private Map<Tag,TaxaDist> tagCntMap=new HashMap<>(200_000_000);

    public FastqToTagCountPlugin2() {
        super(null, false);
    }

    public FastqToTagCountPlugin2(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        countTags(myKeyFile.value(), myEnzyme.value(), myInputDir.value(), myOutputDir.value(), myMaxGoodReads.value(), myMinTagCount.value());
        return null;
    }

    /**
     * Derives a tagCount list for each fastq file in the fastqDirectory.
     *
     * @param keyFileS A key file (a sample key by barcode, with a plate map
     * included).
     * @param enzyme The enzyme used to create the library (currently ApeKI or
     * PstI).
     * @param fastqDirectory Directory containing the fastq files (will be
     * recursively searched).
     * @param outputDir Directory to which the tagCounts files (one per fastq
     * file) will be written.
     * @param maxGoodReads The maximum number of barcoded reads expected in a
     * fastq file
     * @param minCount The minimum number of occurrences of a tag in a fastq
     * file for it to be included in the output tagCounts file
     */
    private void countTags(String keyFileS, String enzyme, String fastqDirectory, String outputDir, int maxGoodReads, int minCount) {
        String[] countFileNames = null;

        File inputDirectory = new File(fastqDirectory);
        System.out.println(inputDirectory.getAbsoluteFile().toString());
        File[] fastqFiles = DirectoryCrawler.listFiles("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$", inputDirectory.getAbsolutePath());
        //                                              (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char' & escape the backslash
        if (fastqFiles.length == 0 || fastqFiles == null) {
            myLogger.warn("Couldn't find any files that end with \".fq\", \".fq.gz\", \".fastq\", \"_fastq.txt\", \"_fastq.gz\", \"_fastq.txt.gz\", \"_sequence.txt\", or \"_sequence.txt.gz\" in the supplied directory.");
            return;
        } else {
            myLogger.info("Using the following FASTQ files:");
            countFileNames = new String[fastqFiles.length];
            for (int i = 0; i < fastqFiles.length; i++) {
                countFileNames[i] = fastqFiles[i].getName().replaceAll("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", ".cnt");
                //                                                                  \\. escape . so it doesn't mean 'any char' & escape the backslash                
                myLogger.info(fastqFiles[i].getAbsolutePath());
            }
        }

        for (int laneNum = 0; laneNum < fastqFiles.length; laneNum++) {
            System.out.println("REading:"+outputDir + File.separator + countFileNames[laneNum]);
            File outputFile = new File(outputDir + File.separator + countFileNames[laneNum]);
            if (outputFile.isFile()) {
                myLogger.warn("An output file " + countFileNames[laneNum] + "\n"
                        + " already exists in the output directory for file " + fastqFiles[laneNum] + ".  Skipping.");
                continue;
            }

            myLogger.info("Reading FASTQ file: " + fastqFiles[laneNum]);
            String[] filenameField = fastqFiles[laneNum].getName().split("_");

            //todo ParseBarcodeRead should be processing barcodes with trie
            ParseBarcodeRead2 thePBR;  // this reads the key file and store the expected barcodes for this lane
            if (filenameField.length == 3) {
                thePBR = new ParseBarcodeRead2(keyFileS, enzyme, filenameField[0], filenameField[1]);
            } else if (filenameField.length == 4) {
                thePBR = new ParseBarcodeRead2(keyFileS, enzyme, filenameField[0], filenameField[2]);
            } // B08AAABXX_s_1_sequence.txt.gz
            else if (filenameField.length == 5) {
                thePBR = new ParseBarcodeRead2(keyFileS, enzyme, filenameField[1], filenameField[3]);
            } else {
                myLogger.error("Error in parsing file name: " + fastqFiles[laneNum]);
                myLogger.error("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
                myLogger.error("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
                continue;
            }
            BiMap<String, Integer> taxaNameToUniqueIndexMap=thePBR.getTaxaNameToUniqueIndexMap();  //This needs to be done once with master file
            int maxTaxaNumber=taxaNameToUniqueIndexMap.size();
            myLogger.info("Total barcodes found in lane:" + thePBR.getBarCodeCount());
            if (thePBR.getBarCodeCount() == 0) {
                myLogger.warn("No barcodes found.  Skipping this flowcell lane.");
                continue;
            }
            String[] taxaNames = new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i] = thePBR.getTheBarcodes(i).getTaxaName();
            }

            int goodBarcodedReads = 0;
            try {
                //todo parallelize this multiple cores
                BufferedReader br = Utils.getBufferedReader(fastqFiles[laneNum], 1<<22);
                int currLine = 0;
                long time=System.nanoTime();
                int allReads = 0;
                goodBarcodedReads = 0;
                String sequence = "";
                String qualityScore = "";
                String temp = br.readLine();
                while ((temp != null) && goodBarcodedReads < maxGoodReads) {
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
                            ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sequence, qualityScore, true, myMinQualScore.value());
                            if (rr != null) {
                                goodBarcodedReads++;
                                Tag tg=TagBuilder.instance(rr.getRead(), rr.getLength());
                                TaxaDist iC=tagCntMap.get(tg);
                                if(iC==null) iC= new TaxaDist(maxTaxaNumber);
                                tagCntMap.put(tg,iC.increment((short) rr.getTaxonIndex()));
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
                //theTC.collapseCounts();
                System.out.println("tagCntMap size"+tagCntMap.size());
                for (Map.Entry<Tag, TaxaDist> entry : tagCntMap.entrySet()) {
                    if(entry.getKey().seqLength()==64 && entry.getValue().totalDepth()>1000) {
                        System.out.println(entry.toString());
                        System.out.println(entry.getValue().taxaDepthMap());
                    }

                }

                //theTC.writeTagCountFile(outputDir + File.separator + countFileNames[laneNum], FilePacking.Byte, minCount);
                myLogger.info("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
                br.close();


            } catch (Exception e) {
                myLogger.error("Good Barcodes Read: " + goodBarcodedReads);
                e.printStackTrace();
            }
            myLogger.info("Finished reading " + (laneNum + 1) + " of " + fastqFiles.length + " sequence files.");
        }
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(FastqToTagCountPlugin2.class);
    // }

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
    public FastqToTagCountPlugin2 inputDirectory(String value) {
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
    public FastqToTagCountPlugin2 keyFile(String value) {
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
    public FastqToTagCountPlugin2 enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

    /**
     * Max good reads per lane
     *
     * @return Max Good Reads
     */
    public Integer maxGoodReads() {
        return myMaxGoodReads.value();
    }

    /**
     * Set Max Good Reads. Max good reads per lane
     *
     * @param value Max Good Reads
     *
     * @return this plugin
     */
    public FastqToTagCountPlugin2 maxGoodReads(Integer value) {
        myMaxGoodReads = new PluginParameter<>(myMaxGoodReads, value);
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
    public FastqToTagCountPlugin2 minTagCount(Integer value) {
        myMinTagCount = new PluginParameter<>(myMinTagCount, value);
        return this;
    }

    /**
     * Output directory to contain .cnt files (one per FASTQ
     * file (one per FASTQ file)
     *
     * @return Output Directory
     */
    public String outputDirectory() {
        return myOutputDir.value();
    }

    /**
     * Set Output Directory. Output directory to contain .cnt
     * files (one per FASTQ file (one per FASTQ file)
     *
     * @param value Output Directory
     *
     * @return this plugin
     */
    public FastqToTagCountPlugin2 outputDirectory(String value) {
        myOutputDir = new PluginParameter<>(myOutputDir, value);
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
    public FastqToTagCountPlugin2 minimumQualityScore(Integer value) {
        myMinQualScore = new PluginParameter<>(myMinQualScore, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Fastq to Tag Count";
    }

    @Override
    public String getToolTipText() {
        return "Fastq to Tag Count";
    }
}
