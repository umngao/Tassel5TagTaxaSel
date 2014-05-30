/*
 * Plugin for the UNEAK pipeline to convert the Tag Pair file into a (fake) TOPM file for reintegration with the normal GBS pipeline
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.tag.AbstractTags;
import net.maizegenetics.dna.tag.UTagPairs;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.ArgsEngine;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 *
 * @author Jason Wallace
 *
 * This plugin takes the TagPairs file produced by earlier steps in the UNEAK
 * pipeline and converts it to a TOPM alignment file that can be used by
 * DiscoverySNPCallerPlugin to call SNPs as if from a reference genome.
 *
 */
public class UTagPairToTOPMPlugin extends AbstractPlugin {

    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UTagPairToTOPMPlugin.class);

    //Internal data for converting tag pairs
    private UTagPairs tp;   //Structure to read in tag pairs
    byte myStrand = 1, myMultimaps = 1, myDcoP = Byte.MIN_VALUE, myMapP = Byte.MIN_VALUE;   //Dummy values to feed to the TOPM structure
    //int myMaxMapping = 1, myMaxVariants = 8;    //HDF5-specific dummy values for output; not currently in use since DiscoverySNPCallerPlugin doesn't take HDF5-TOPM as input

    //Parameters for this plugin
    private PluginParameter<Integer> chrom
            = new PluginParameter.Builder<>("chrom", 1, Integer.class)
            .description("Chromosome to start numbering at")
            .guiName("Start chromosome")
            .build();
    private PluginParameter<Integer> distance
            = new PluginParameter.Builder<>("distance", 1000, Integer.class)
            .description("Distance to pad between each tag pair")
            .guiName("Pad distance")
            .units("base pairs")
            .build();
    private PluginParameter<String> infile
            = new PluginParameter.Builder<>("input", null, String.class)
            .description("Input file of matched tag pairs")
            .required(true)
            .inFile()
            .build();
    private PluginParameter<String> textOutputFile
            = new PluginParameter.Builder<>("toText", null, String.class)
            .description("File to output TOPM in text format")
            .required(false)
            .outFile()
            .build();
    private PluginParameter<String> binaryOutputFile
            = new PluginParameter.Builder<>("toBinary", null, String.class)
            .description("File to output TOPM in binary format")
            .required(false)
            .outFile()
            .build();


    public UTagPairToTOPMPlugin() {
        super(null, false);
    }

    public UTagPairToTOPMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        //Make TOPM
        tp = new UTagPairs(inputFile());
        TagsOnPhysicalMap myTopm = makeTopmFromTagPairs(tp);

        //Output to specified file(s)
        if (textOutputFile() != null) {
            logger.info("Outputting TOPM in text format to " + textOutputFile());
            myTopm.writeTextFile(new File(textOutputFile()));
        }
        if (binaryOutputFile() != null) {
            logger.info("Outputting TOPM in binary format to " + binaryOutputFile());
            myTopm.writeBinaryFile(new File(binaryOutputFile()));
        }

        return null;
    }

    /**
     * Convert the TagPair object to a TOPM object in memory, filling in things
     * with dummy data where needed
     */
    private TagsOnPhysicalMap makeTopmFromTagPairs(UTagPairs tp) {
        HelperTags tempTags = new HelperTags(tp.getTagNum(), tp.getTag(0).length);
        //Load in all the data for each tag pair; having trouble finding the function to load actual sequences
        for (int i = 0; i < tp.getTagNum(); i++) {
            tempTags.setTag(i, tp.getTag(i), tp.getTagLength(i));
        }
        TagsOnPhysicalMap myTopm = new TagsOnPhysicalMap(tempTags);
        int currPos = 1;
        int currChrom = startChrom();
        for (int i = 0; i < tp.getTagNum(); i++) {

            myTopm.setChromoPosition(i, currChrom, myStrand, currPos, currPos + tp.getTagLength(i) - 1);
            myTopm.setDivergence(i, (byte) 0);

            //These may not be necessary; don't know
            myTopm.setDcoP(i, Byte.MIN_VALUE);  //May have to alter TOPM class to do this
            myTopm.setMapP(i, Byte.MIN_VALUE);
            myTopm.setMultimaps(i, (byte) 1);

            //Increment position after odd-numbered tags (so pairs are at the same position)
            if (i % 2 == 1) {
                currPos += padDistance();
            }
            //If over max interger value, increment chromosome and start over
            if (currPos >= Integer.MAX_VALUE) {
                currChrom++;
                currPos = 1;
            }
        }
        return myTopm;
    }

    @Override
    protected void postProcessParameters(){
        //Test that at least one output file supplied
        if((binaryOutputFile() == null) && (textOutputFile() == null)){
            throw new IllegalArgumentException("\n\nMust specify at least one output file (text or binary).\n\n");
        }
    }

    @Override
    public String pluginDescription(){
        return "This plugin takes a tag-pair file (from UTagCountToTagPairPlugin) and converts it into a TOPM (tags on physical map) " +
                "file for use in the GBS pipeline. The resulting chromosome coordinates and other data are just filler to comply with " +
                "the TOPM file specifications.";
    }

    //Parameter get/set functions
    public UTagPairToTOPMPlugin inputFile(String filename){
        setParameter(infile.cmdLineName(), filename);
        return this;
    }

    public UTagPairToTOPMPlugin textOutfile(String filename){
        setParameter(textOutputFile.cmdLineName(), filename);
        return this;
    }

    public UTagPairToTOPMPlugin binaryOutfile(String filename){
        setParameter(binaryOutputFile.cmdLineName(), filename);
        return this;
    }

    public UTagPairToTOPMPlugin startChrom(Integer value){
        setParameter(chrom.cmdLineName(), value);
        return this;
    }

    public UTagPairToTOPMPlugin padDistance(Integer value){
        setParameter(distance.cmdLineName(), value);
        return this;
    }

    public String inputFile(){
        return infile.value();
    }

    public String textOutputFile(){
        return textOutputFile.value();
    }

    public String binaryOutputFile(){
        return binaryOutputFile.value();
    }

    public Integer startChrom(){
        return chrom.value();
    }

    public Integer padDistance(){
        return distance.value();
    }

    @Override
    public ImageIcon getIcon() {
        return(null);
    }

    @Override
    public String getButtonName() {
        return("Tag Pairs to TOPM");
    }

    @Override
    public String getToolTipText() {
        return("Tag Pairs to TOPM");
    }
}

/*
 * A small helper class that exists solely to pass tag info to TOPM
 */
class HelperTags extends AbstractTags {

    /* Inherited
     protected int tagLengthInLong;  //
     protected long[][] tags;  //Index in rows, actual tag components in columns
     protected byte[] tagLength;  // length of tag (number of bases)  // 1 byte
     */
    public HelperTags(int numtags, int myTagLengthInLong) {
        tagLengthInLong = myTagLengthInLong;
        tags = new long[tagLengthInLong][numtags];
        tagLength = new byte[numtags];
    }

    public void setTag(int index, long[] tagValue, byte myTagLength) {
        tagLength[index] = myTagLength;
        for (int i = 0; i < tagValue.length; i++) {
            tags[i][index] = tagValue[i];
        }
    }
}
