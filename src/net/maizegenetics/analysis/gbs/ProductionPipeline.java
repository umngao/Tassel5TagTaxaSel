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
import java.io.FilenameFilter;

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
    private static final String READY_FILE_NAME = "ready.txt";
    private static final String LOCK_FILE_NAME = "lock.txt";

    public enum PARAMETERS {

        inputDirectory, keyFile, enzyme, productionTOPM, outputGenotypeFile, archiveDirectory
    };

    protected PluginParameter<String> myInputDirectory = new PluginParameter.Builder<>(PARAMETERS.inputDirectory, null, String.class).required(true).inDir()
            .description("Input directory containing subdirectories with fastq AND/OR qseq files").build();
    protected PluginParameter<String> myEnzyme = new PluginParameter.Builder<>(PARAMETERS.enzyme, null, String.class).required(true)
            .description("Enzyme used to create the GBS library").build();
    protected PluginParameter<String> myProductionTOPM = new PluginParameter.Builder<>(PARAMETERS.productionTOPM, null, String.class).required(true).inFile()
            .description("Physical map file containing tags and corresponding variants (production TOPM)").build();
    protected PluginParameter<String> myOutputGenotypeFile = new PluginParameter.Builder<>(PARAMETERS.outputGenotypeFile, null, String.class).required(true).outFile()
            .description("Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)").build();
    protected PluginParameter<String> myArchiveDirectory = new PluginParameter.Builder<>(PARAMETERS.archiveDirectory, null, String.class).required(true).outDir()
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

            String lockFilename = inputDirectory() + File.separator + LOCK_FILE_NAME;
            if (new File(lockFilename).exists()) {
                myLogger.warn("Production Pipeline already running.  File exists: " + lockFilename + "  Aborting...");
                return null;
            }

            File inputDirectory = new File(inputDirectory());
            String[] directories = inputDirectory.list(new FilenameFilter() {
                @Override
                public boolean accept(File current, String name) {
                    return new File(current, name).isDirectory();
                }
            });

            for (String current : directories) {
                String fullDirName = inputDirectory() + File.separator + current;
                processSubDirectory(fullDirName);
                File currentLocation = new File(fullDirName);
                File newLocation = new File(archiveDirectory() + current);
                currentLocation.renameTo(newLocation);
                myLogger.info("Moved : " + currentLocation.getAbsolutePath() + " to: " + newLocation.getAbsolutePath());
            }

            return null;
        } finally {
            LoggingUtils.closeLogfile();
        }

    }

    private void processSubDirectory(String subDirectory) {

        String readyFilename = subDirectory + File.separator + READY_FILE_NAME;
        File readyFile = new File(readyFilename);
        if (readyFile.exists()) {
            myLogger.info("Processing directory: " + subDirectory);
            readyFile.delete();
        } else {
            myLogger.warn("This directory is not ready yet: " + subDirectory);
            return;
        }

        String keyFile = subDirectory + File.separator + Utils.getFilename(subDirectory) + ".key";
        if (!new File(keyFile).exists()) {
            myLogger.error("Keyfile doesn't exist: " + keyFile);
            return;
        }

        String[] rawSeqFileNames = DirectoryCrawler.listFileNames(ProductionSNPCallerPlugin.rawSeqFileNameRegex, subDirectory);
        if ((rawSeqFileNames == null) || (rawSeqFileNames.length == 0)) {
            myLogger.warn("No sequence files in directory: " + subDirectory);
            return;
        }

        myLogger.info("Raw Sequence Files: " + Arrays.deepToString(rawSeqFileNames));
        myLogger.info("Parameters Passed to ProductionSNPCallerPlugin: " + Arrays.deepToString(getPluginArgs(subDirectory, keyFile)));

        ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin();
        plugin.setParameters(getPluginArgs(subDirectory, keyFile));

        printParameterValues();
        plugin.performFunction(null);

    }

    private String[] getPluginArgs(String inputDir, String keyFile) {
        String[] args = {
            "-i", inputDir,
            "-k", keyFile,
            "-e", enzyme(),
            "-o", outputGenotypeFile(),
            "-m", productionTOPM()
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
     * Input directory containing subdirectories with fastq AND/OR qseq files
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDirectory.value();
    }

    /**
     * Set Input Directory. Input directory containing subdirectories with fastq
     * AND/OR qseq files
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
