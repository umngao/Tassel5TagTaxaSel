/*
 * BinaryToTextPlugin
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.collect.Range;

import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.Arrays;

/**
 *
 * @author Terry Casstevens
 */
public class BinaryToTextPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(BinaryToTextPlugin.class);

    public static enum FILE_TYPES {

        TOPM, TagCounts, TBTByte

    };

    private PluginParameter<String> myInputFile = new PluginParameter.Builder<String>("i", null, String.class).guiName("Input File").required(true).inFile().build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile().build();
    private PluginParameter<FILE_TYPES> myFileType = new PluginParameter.Builder<FILE_TYPES>("t", FILE_TYPES.TOPM, FILE_TYPES.class).guiName("File Type").range(Range.encloseAll(Arrays.asList(FILE_TYPES.values()))).build();

    public BinaryToTextPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        switch (fileType()) {
            case TOPM:
                TagsOnPhysicalMap topm = new TagsOnPhysicalMap(inputFile(), true);
                topm.writeTextFile(new File(outputFile()));
                break;
            case TagCounts:
                TagCounts tc = new TagCounts(inputFile(), FilePacking.Byte);
                tc.writeTagCountFile(outputFile(), FilePacking.Text, 0);
                break;
            case TBTByte:
                TagsByTaxaByte tbtbyte = new TagsByTaxaByte(inputFile(), FilePacking.Byte);
                tbtbyte.writeDistFile(new File(outputFile()), FilePacking.Text, 0);
                break;
        }

        return null;

    }

    public BinaryToTextPlugin inputFile(String filename) {
        setParameter(myInputFile.cmdLineName(), filename);
        return this;
    }

    public String inputFile() {
        return myInputFile.value();
    }

    public BinaryToTextPlugin outputFile(String filename) {
        setParameter(myOutputFile.cmdLineName(), filename);
        return this;
    }

    public String outputFile() {
        return myOutputFile.value();
    }

    public BinaryToTextPlugin fileType(FILE_TYPES type) {
        setParameter(myFileType.cmdLineName(), type.toString());
        return this;
    }

    public FILE_TYPES fileType() {
        return myFileType.value();
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Binary to Text";
    }

    @Override
    public String getToolTipText() {
        return "Binary to Text";
    }
}
