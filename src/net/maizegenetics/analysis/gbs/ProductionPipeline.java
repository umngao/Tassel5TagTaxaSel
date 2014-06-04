/*
 * ProductionPipeline
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.util.LoggingUtils;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;

import java.io.File;

import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

/**
 *
 * @author Terry Casstevens
 */
public class ProductionPipeline extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ProductionPipeline.class);
    private static final SimpleDateFormat LOGGING_DATE_FORMAT = new SimpleDateFormat("yyyyMMdd HH:mm:ss");

    public enum PARAMETERS {

        inputDirectory, keyFile, enzyme, productionTOPM, outputGenotypeFile, archiveDirectory
    };

    protected PluginParameter<String> myInputDirectory = new PluginParameter.Builder<String>(PARAMETERS.inputDirectory, null, String.class).required(true).inDir()
            .description("Input directory containing fastq AND/OR qseq files").build();
    protected PluginParameter<String> myKeyFile = new PluginParameter.Builder<String>(PARAMETERS.keyFile, null, String.class).required(true).inFile()
            .description("Barcode Key File").build();
    protected PluginParameter<String> myEnzyme = new PluginParameter.Builder<String>(PARAMETERS.enzyme, null, String.class).required(true)
            .description("Enzyme used to create the GBS library").build();
    protected PluginParameter<String> myProductionTOPM = new PluginParameter.Builder<String>(PARAMETERS.productionTOPM, null, String.class).required(true).inFile()
            .description("Physical map file containing tags and corresponding variants (production TOPM)").build();
    protected PluginParameter<String> myOutputGenotypeFile = new PluginParameter.Builder<String>(PARAMETERS.outputGenotypeFile, null, String.class).required(true).outFile()
            .description("Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)").build();
    protected PluginParameter<String> myArchiveDirectory = new PluginParameter.Builder<String>(PARAMETERS.archiveDirectory, null, String.class).required(true).outDir()
            .description("Archive directory where to move processed files").build();

    private String myOutputDirectory;

    public ProductionPipeline(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        if (myOutputGenotypeFile.value() != null) {
            myOutputDirectory = Utils.getDirectory(myOutputGenotypeFile.value());
            setupLogfile();
        }
        myLogger.info(getTimeStamp());
    }

    @Override
    public DataSet processData(DataSet input) {

        try {
            String[] rawSeqFileNames = DirectoryCrawler.listFileNames(ProductionSNPCallerPlugin.rawSeqFileNameRegex, myInputDirectory.value());
            if ((rawSeqFileNames == null) || (rawSeqFileNames.length == 0)) {
                return null;
            }

            myLogger.info("Raw Sequence Files: " + Arrays.deepToString(rawSeqFileNames));
            myLogger.info("Parameters Passed to ProductionSNPCallerPlugin: " + Arrays.deepToString(getPluginArgs()));

            ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin();
            plugin.setParameters(getPluginArgs());

            printParameterValues();
            plugin.performFunction(null);

            for (int i = 0; i < rawSeqFileNames.length; i++) {
                File currentLocation = new File(rawSeqFileNames[i]);
                File newLocation = new File(myArchiveDirectory.value() + Utils.getFilename(rawSeqFileNames[i]));
                currentLocation.renameTo(newLocation);
                myLogger.info("Moved : " + currentLocation.getAbsolutePath() + " to: " + newLocation.getAbsolutePath());
            }

            return null;
        } finally {
            LoggingUtils.closeLogfile();
        }

    }

    private String[] getPluginArgs() {
        String[] args = {
            "-i", myInputDirectory.value(),
            "-k", myKeyFile.value(),
            "-e", myEnzyme.value(),
            "-o", myOutputGenotypeFile.value(),
            "-m", myProductionTOPM.value()
        };
        return args;
    }

    private void setupLogfile() {

        String todayDate = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String logFileName = todayDate + "_" + "ProductionPipeline" + ".log";
        logFileName = myOutputDirectory + "/" + logFileName;

        myLogger.info("Log File: " + logFileName);
        try {
            LoggingUtils.setupLogfile(logFileName);
        } catch (Exception e) {
            throw new IllegalArgumentException("ProductionPipeline: setupLogfile: " + logFileName + " doesn't exist.");
        }

    }

    /**
     * Convenience method to provide uniformly labeled timestamps
     */
    private static String getTimeStamp() {
        return "Timestamp: " + LOGGING_DATE_FORMAT.format(new Date()) + ": ";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Production Pipeline";
    }

    @Override
    public String getToolTipText() {
        return "Production Pipeline";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ProductionPipeline.class);
    // }
    /**
     * Input directory containing fastq AND/OR qseq files
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDirectory.value();
    }

    /**
     * Set Input Directory. Input directory containing fastq AND/OR qseq files
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public ProductionPipeline inputDirectory(String value) {
        myInputDirectory = new PluginParameter<>(myInputDirectory, value);
        return this;
    }

    /**
     * Barcode Key File
     *
     * @return Key File
     */
    public String keyFile() {
        return myKeyFile.value();
    }

    /**
     * Set Key File. Barcode Key File
     *
     * @param value Key File
     *
     * @return this plugin
     */
    public ProductionPipeline keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }

    /**
     * Enzyme used to create the GBS library
     *
     * @return Enzyme
     */
    public String enzyme() {
        return myEnzyme.value();
    }

    /**
     * Set Enzyme. Enzyme used to create the GBS library
     *
     * @param value Enzyme
     *
     * @return this plugin
     */
    public ProductionPipeline enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

    /**
     * Physical map file containing tags and corresponding variants (production
     * TOPM)
     *
     * @return Production T O P M
     */
    public String productionTOPM() {
        return myProductionTOPM.value();
    }

    /**
     * Set Production T O P M. Physical map file containing tags and
     * corresponding variants (production TOPM)
     *
     * @param value Production T O P M
     *
     * @return this plugin
     */
    public ProductionPipeline productionTOPM(String value) {
        myProductionTOPM = new PluginParameter<>(myProductionTOPM, value);
        return this;
    }

    /**
     * Output (target) HDF5 genotypes file to add new genotypes to (new file
     * created if it doesn't exist)
     *
     * @return Output Genotype File
     */
    public String outputGenotypeFile() {
        return myOutputGenotypeFile.value();
    }

    /**
     * Set Output Genotype File. Output (target) HDF5 genotypes file to add new
     * genotypes to (new file created if it doesn't exist)
     *
     * @param value Output Genotype File
     *
     * @return this plugin
     */
    public ProductionPipeline outputGenotypeFile(String value) {
        myOutputGenotypeFile = new PluginParameter<>(myOutputGenotypeFile, value);
        return this;
    }

    /**
     * Archive directory where to move processed files
     *
     * @return Archive Directory
     */
    public String archiveDirectory() {
        return myArchiveDirectory.value();
    }

    /**
     * Set Archive Directory. Archive directory where to move processed files
     *
     * @param value Archive Directory
     *
     * @return this plugin
     */
    public ProductionPipeline archiveDirectory(String value) {
        myArchiveDirectory = new PluginParameter<>(myArchiveDirectory, value);
        return this;
    }
}
