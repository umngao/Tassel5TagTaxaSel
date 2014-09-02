/*
 * SAMConverterPlugin
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;

import java.io.*;
import java.util.Set;

import javax.swing.ImageIcon;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 * This class can read in a CBSU TagMapFile into the gbs.TagsOnPhysicalMap data
 * structure.
 *
 * @author harriman
 *
 */
public final class SAMToGBSdbPlugin extends AbstractPlugin {

    boolean cleanCutSites = true;
    private static final Logger myLogger = Logger.getLogger(SAMToGBSdbPlugin.class);

    private PluginParameter<String> myInputFile = new PluginParameter.Builder<String>("i", null, String.class).guiName("SAM Input File").required(true).inFile()
            .description("Name of input file in SAM text format").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("GBS DB File").required(true).outFile()
            .description("Name of output file (Default: output.topm.bin)").build();

    public SAMToGBSdbPlugin() {
        super(null, false);
    }

    public SAMToGBSdbPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        final int name = 0, flag = 1, chr = 2, pos = 3, cigar = 5, tagS = 9, alignScoreIndex=11; // column indices in inputLine
        try {
            BufferedReader bw=Utils.getBufferedReader(sAMInputFile());
            TagDataWriter tagData=new TagDataSQLite(gBSDBFile());
            Set<Tag> knownTags=tagData.getTags();
            Multimap<Tag,Position> tagPositions= HashMultimap.create(knownTags.size(),2);
            String inputLine;
            while((inputLine=bw.readLine())!=null) {
                if(inputLine.startsWith("@")) continue;
                String[] s=inputLine.split("\\s");
                Tag tag= TagBuilder.instance(s[tagS]);
                Chromosome chromosome=new Chromosome(s[chr].replace("chr", ""));
                String alignmentScore=s[alignScoreIndex].split(":")[2];
                Position position=new GeneralPosition
                        .Builder(chromosome,Integer.parseInt(s[pos]))
                        .addAnno("mappingapproach", "Bowtie")
                        .addAnno("cigar", s[cigar])
                        .addAnno("supportvalue", alignmentScore)  //todo include again
                        .build();
                //System.out.println(inputLine);
                tagPositions.put(tag,position);
            }
            bw.close();
            tagData.putTagAlignments(tagPositions);

//            myLogger.info("Finished converting binary tag count file to fastq."
//                    + "\nTotal number of tags written: " + count.get() + " (above minCount of " + minCount() + ")"
//                    + "\nOuput fastq file: " + outputFile() + "\n\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        return null;
    }

    /**
     * Reads SAM files output from BWA or bowtie2
     */
//    private void readSAMFile(String inputFileName, int tagLengthInLong) {
//        System.out.println("Reading SAM format tag alignment from: " + inputFileName);
//        this.tagLengthInLong = tagLengthInLong;
//        String inputStr = "Nothing has been read from the file yet";
//        int nHeaderLines = countTagsInSAMfile(inputFileName); // detects if the file is bowtie2, initializes topm matrices
//        int tagIndex = Integer.MIN_VALUE;
//        try {
//            BufferedReader br;
//            if (inputFileName.endsWith(".gz")) {
//                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
//            } else {
//                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
//            }
//            for (int i = 0; i < nHeaderLines; i++) {
//                br.readLine();
//            } // Skip over the header
//            for (tagIndex = 0; tagIndex < myNumTags; tagIndex++) {
//                inputStr = br.readLine();
//                parseSAMAlignment(inputStr, tagIndex);
//                if (tagIndex % 1000000 == 0) {
//                    System.out.println("Read " + tagIndex + " tags.");
//                }
//            }
//            br.close();
//        } catch (Exception e) {
//            System.out.println("\n\nCatch in reading SAM alignment file at tag " + tagIndex + ":\n\t" + inputStr + "\nError: " + e + "\n\n");
//            e.printStackTrace();
//            System.exit(1);
//        }
//    }
//
//    private int countTagsInSAMfile(String inputFileName) {
//        mySAMFormat = SAMFormat.BWA;  // format is BWA by default
//        myNumTags = 0;
//        int nHeaderLines = 0;
//        String currLine = null;
//        try {
//            String[] inputLine;
//            ArrayList<String> chrNames = new ArrayList<String>();
//            BufferedReader br;
//            if (inputFileName.endsWith(".gz")) {
//                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
//            } else {
//                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
//            }
//            while ((currLine = br.readLine()) != null) {
//                inputLine = currLine.split("\\s");
//                if (inputLine[0].contains("@")) {
//                    //SAM files produced by Bowtie2 contain the string "@PG     ID:bowtie2      PN:bowtie2 "
//                    if (inputLine[1].contains("bowtie2")) {
//                        mySAMFormat = SAMFormat.BOWTIE2;
//                    }
//                    nHeaderLines++;
//                } else {
//                    String chr = inputLine[2];
//                    if (!chrNames.contains(chr)) {
//                        chrNames.add(chr);
//                    }
//                    myNumTags++;
//                    if (myNumTags % 1000000 == 0) {
//                        System.out.println("Counted " + myNumTags + " tags.");
//                    }
//                }
//            }
//            br.close();
//            System.out.println("Found " + myNumTags + " tags in SAM file.  Assuming " + mySAMFormat + " file format.");
//        } catch (Exception e) {
//            System.out.println("Catch in counting lines of alignment file at line " + currLine + ": " + e);
//            e.printStackTrace();
//            System.exit(1);
//        }
//        initMatrices(myNumTags);
//        return nHeaderLines;
//    }
//
//
//    private void parseSAMAlignment(String inputStr, int tagIndex) {
//        String[] inputLine = inputStr.split("\t");
//        int name = 0, flag = 1, chr = 2, pos = 3, cigar = 5, tagS = 9; // column indices in inputLine
//        String nullS = this.getNullTag();
//        byte currStrand = ((Integer.parseInt(inputLine[flag]) & 16) == 16) ? (byte) -1 : (byte) 1; // bit 0x10 (= 2^4 = 16) is set: REVERSE COMPLEMENTED
//        if ((Integer.parseInt(inputLine[flag]) & 4) == 4) {  // bit 0x4 (= 2^2 = 4) is set: NO ALIGNMENT
//            recordLackOfSAMAlign(tagIndex, inputLine[tagS], inputLine[name], nullS, currStrand);
//        } else {  // aligns to one or more positions
//            HashMap<String, Integer> SAMFields = parseOptionalFieldsFromSAMAlignment(inputLine);
//            byte bestHits = (byte) Math.min(SAMFields.get("nBestHits"), Byte.MAX_VALUE);
//            byte editDist = (byte) Math.min(SAMFields.get("editDist"), Byte.MAX_VALUE);
//            recordSAMAlign(
//                    tagIndex,
//                    inputLine[tagS],
//                    inputLine[name],
//                    nullS,
//                    bestHits,
//                    inputLine[chr],
//                    currStrand,
//                    Integer.parseInt(inputLine[pos]),
//                    inputLine[cigar],
//                    editDist);
//        }
//    }
//
//    private HashMap<String, Integer> parseOptionalFieldsFromSAMAlignment(String[] inputLine) {
//        HashMap<String, Integer> SAMFields = new HashMap<String, Integer>();
//        if (mySAMFormat == SAMFormat.BWA) {
//            for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
//                if (inputLine[field].regionMatches(0, "X0", 0, 2)) {        // X0 = SAM format for # of "high-quality" alignments of this query.  Specific to BWA.
//                    SAMFields.put("nBestHits", Integer.parseInt(inputLine[field].split(":")[2]));
//                } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
//                    SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
//                }
//            }
//        } else {  // bowtie2 -M format
//            for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
//                if (inputLine[field].regionMatches(0, "AS", 0, 2)) {        // AS = SAM format for alignment score of the best alignment.  Specific to bowtie2.
//                    SAMFields.put("bestScore", Integer.parseInt(inputLine[field].split(":")[2]));
//                } else if (inputLine[field].regionMatches(0, "XS", 0, 2)) { // XS = SAM format for alignment score of 2nd best alignment.  Specific to bowtie2.
//                    SAMFields.put("nextScore", Integer.parseInt(inputLine[field].split(":")[2]));
//                } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
//                    SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
//                }
//            }
//            if (SAMFields.containsKey("bestScore")) {
//                if (SAMFields.containsKey("nextScore")) {
//                    if (SAMFields.get("bestScore") > SAMFields.get("nextScore")) {
//                        SAMFields.put("nBestHits", 1);
//                    } else {
//                        SAMFields.put("nBestHits", 99);  // 99 will stand for an unknown # of multiple hits
//                    }
//                } else {
//                    SAMFields.put("nBestHits", 1);
//                }
//            }
//        }
//        return SAMFields;
//    }
//
//    private void writeLogFile(TagsOnPhysicalMap topm) {
//        try {
//            DataOutputStream report = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile() + ".log"), 65536));
//            int[] aligned = topm.mappedTags();
//            int unique = 0, multi = 1;  // the indices of aligned
//            int unaligned = topm.getTagCount() - aligned[unique] - aligned[multi];
//            report.writeBytes(
//                    "Input file: " + inputFile() + "\n"
//                    + "Output file: " + outputFile() + "\n"
//                    + "Total " + topm.getTagCount() + " tags\n\t"
//                    + aligned[unique] + " were aligned to unique postions\n\t"
//                    + aligned[multi] + " were aligned to multiple postions\n\t"
//                    + unaligned + " could not be aligned.\n\n");
//            int[] dist = topm.mappingDistribution();
//            report.writeBytes("nPositions  nTags\n");
//            for (int i = 0; i < dist.length; i++) {
//                if (dist[i] > 0) {
//                    if (i < 10) {
//                        report.writeBytes(i + "           " + dist[i] + "\n");
//                    } else if (i < 100) {
//                        report.writeBytes(i + "          " + dist[i] + "\n");
//                    } else if (i < 1000) {
//                        report.writeBytes(i + "         " + dist[i] + "\n");
//                    }
//                }
//            }
//            report.close();
//        } catch (Exception e) {
//            myLogger.warn("Caught exception while writing log file: " + e);
//        }
//    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(SAMToGBSdbPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TagData runPlugin(DataSet input) {
        return (TagData) performFunction(input).getData(0).getData();
    }

    /**
     * Name of input file in SAM text format
     *
     * @return SAM Input File
     */
    public String sAMInputFile() {
        return myInputFile.value();
    }

    /**
     * Set SAM Input File. Name of input file in SAM text
     * format
     *
     * @param value SAM Input File
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin sAMInputFile(String value) {
        myInputFile = new PluginParameter<>(myInputFile, value);
        return this;
    }

    /**
     * Name of output file (Default: output.topm.bin)
     *
     * @return GBS DB File
     */
    public String gBSDBFile() {
        return myOutputFile.value();
    }

    /**
     * Set GBS DB File. Name of output file (Default: output.topm.bin)
     *
     * @param value GBS DB File
     *
     * @return this plugin
     */
    public SAMToGBSdbPlugin gBSDBFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }


    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "SAM to TOPM Converter";
    }

    @Override
    public String getToolTipText() {
        return "SAM to TOPM Converter";
    }
}
