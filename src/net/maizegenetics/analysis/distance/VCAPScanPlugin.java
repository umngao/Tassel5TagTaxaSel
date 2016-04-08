/*
 *  VCAPScanPlugin
 * 
 *  Created on Apr 4, 2016
 */
package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.DefaultPluginListener;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.plugindef.ThreadedPluginListener;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class VCAPScanPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(VCAPScanPlugin.class);

    public static enum SCAN_METHOD {

        Chromosome,
        Site_Blocks,
        Bed_File
    };

    private PluginParameter<SCAN_METHOD> myMethod = new PluginParameter.Builder<>("method", SCAN_METHOD.Chromosome, SCAN_METHOD.class)
            .guiName("Scan method")
            .range(SCAN_METHOD.values())
            .description("")
            .build();

    private PluginParameter<Integer> myNumSitesPerBlock = new PluginParameter.Builder<>("numSitesPerBlock", 10000, Integer.class)
            .description("")
            .build();

    private PluginParameter<String> myBedFile = new PluginParameter.Builder<>("bedFile", null, String.class)
            .description("")
            .inFile()
            .build();

    private PluginParameter<String> myPhenotypeFile = new PluginParameter.Builder<>("phenotypeFile", null, String.class)
            .description("")
            .inFile()
            .required(true)
            .build();

    private PluginParameter<String> myLDAKCommand = new PluginParameter.Builder<>("ldakCommand", "ldak.4.9.fast", String.class)
            .description("")
            .build();

    private PluginParameter<String> myOutputDir = new PluginParameter.Builder<>("outputDir", ".", String.class)
            .description("")
            .outDir()
            .build();

    private PluginParameter<String> myWholeMatrix = new PluginParameter.Builder<>("wholeMatrix", null, String.class)
            .description("")
            .inFile()
            .build();

    public VCAPScanPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if ((alignInList == null) || (alignInList.size() != 1)) {
            throw new IllegalArgumentException("VCAPScanPlugin: Must input exactly one genotype table.");
        }
    }

    @Override
    protected void postProcessParameters() {
        if (outputDir().charAt(outputDir().length() - 1) != '/') {
            outputDir(outputDir() + "/");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        GenotypeTable genotype = (GenotypeTable) input.getDataOfType(GenotypeTable.class).get(0).getData();
        List<String> matrixFiles = new ArrayList<>();

        switch (method()) {

            case Chromosome:
                ForkJoinPool threadPool = new ForkJoinPool();

                KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                        .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
                kinshipPlugin.addListener(DefaultPluginListener.getInstance());

                SubtractDistanceMatrixPlugin subtractPlugin = new SubtractDistanceMatrixPlugin(null, false)
                        .wholeMatrix(createWholeMatrixIfNeeded(input));

                for (Chromosome currentChr : genotype.chromosomes()) {

                    String chrStr = currentChr.getName();
                    FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                            .startChr(currentChr)
                            .endChr(currentChr);
                    DataSet genotypeDataSet = filterChr.performFunction(input);

                    DataSet part = kinshipPlugin.performFunction(genotypeDataSet);

                    DataSet rest = subtractPlugin.performFunction(part);

                    ExportPlugin exportPlugin = new ExportPlugin(null, false);
                    exportPlugin.setAlignmentFileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                    String saveFilename = outputDir() + "Kinship_" + chrStr;
                    matrixFiles.add(saveFilename);
                    exportPlugin.setSaveFile(saveFilename);
                    threadPool.submit(new ThreadedPluginListener(exportPlugin, new PluginEvent(new DataSet(part.getData(0), part.getCreator()))));

                    ExportPlugin exportPlugin1 = new ExportPlugin(null, false);
                    exportPlugin1.setAlignmentFileType(FileLoadPlugin.TasselFileType.SqrMatrixBin);
                    exportPlugin1.setSaveFile(saveFilename + "Rest");
                    threadPool.submit(new ThreadedPluginListener(exportPlugin1, new PluginEvent(new DataSet(rest.getData(0), part.getCreator()))));

                }

                threadPool.shutdown();
                try {
                    threadPool.awaitTermination(20, TimeUnit.MINUTES);
                } catch (Exception e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException("VCAPScanPlugin: processData: problem: " + e.getMessage());
                }
                break;

            case Site_Blocks:
                throw new UnsupportedOperationException();

            case Bed_File:
                throw new UnsupportedOperationException();

        }

        runLDAK(matrixFiles);

        return null;

    }

    private String createWholeMatrixIfNeeded(DataSet input) {
        if ((wholeMatrix() == null) || (wholeMatrix().isEmpty())) {
            KinshipPlugin kinshipPlugin = new KinshipPlugin(null, false)
                    .kinshipMethod(KinshipPlugin.KINSHIP_METHOD.Centered_IBS);
            kinshipPlugin.addListener(DefaultPluginListener.getInstance());
            DataSet whole = kinshipPlugin.performFunction(input);

            ExportPlugin exportPlugin = new ExportPlugin(null, false);
            exportPlugin.setAlignmentFileType(FileLoadPlugin.TasselFileType.SqrMatrix);
            String saveFilename = outputDir() + "kinship_whole.txt";
            exportPlugin.setSaveFile(saveFilename);
            exportPlugin.performFunction(whole);
            return saveFilename;
        } else {
            return wholeMatrix();
        }
    }

    private void runLDAK(List<String> matrixFiles) {
        
        // ldak.4.9.fast --reml results_chr1 --mgrm kinship_list.txt --pheno NAM_ap_dta_multiblup.txt --kinship-details NO
        for (String filename : matrixFiles) {
            
            try (BufferedWriter writer = Utils.getBufferedWriter("kinship_list.txt")) {
                writer.write(filename);
                writer.write("\n");
                writer.write(filename + "Rest");
                writer.write("\n");
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
            }

            String command = ldakCommand()
                    + " --reml " + "results" + Utils.getFilename(filename)
                    + " --mgrm " + "kinship_list.txt"
                    + " --pheno " + phenotypeFile()
                    + " --kinship-details NO";
            
            myLogger.info("command: " + command);

            try {
                Process process = Runtime.getRuntime().exec(command);
                BufferedReader reader = new BufferedReader(
                        new InputStreamReader(process.getInputStream()));
                String line = null;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
                process.waitFor();
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
            }

        }

    }

    /**
     * Scan method
     *
     * @return Scan method
     */
    public SCAN_METHOD method() {
        return myMethod.value();
    }

    /**
     * Set Scan method. Scan method
     *
     * @param value Scan method
     *
     * @return this plugin
     */
    public VCAPScanPlugin method(SCAN_METHOD value) {
        myMethod = new PluginParameter<>(myMethod, value);
        return this;
    }

    /**
     * Num Sites Per Block
     *
     * @return Num Sites Per Block
     */
    public Integer numSitesPerBlock() {
        return myNumSitesPerBlock.value();
    }

    /**
     * Set Num Sites Per Block. Num Sites Per Block
     *
     * @param value Num Sites Per Block
     *
     * @return this plugin
     */
    public VCAPScanPlugin numSitesPerBlock(Integer value) {
        myNumSitesPerBlock = new PluginParameter<>(myNumSitesPerBlock, value);
        return this;
    }

    /**
     * Bed File
     *
     * @return Bed File
     */
    public String bedFile() {
        return myBedFile.value();
    }

    /**
     * Set Bed File. Bed File
     *
     * @param value Bed File
     *
     * @return this plugin
     */
    public VCAPScanPlugin bedFile(String value) {
        myBedFile = new PluginParameter<>(myBedFile, value);
        return this;
    }

    /**
     * Phenotype File
     *
     * @return Phenotype File
     */
    public String phenotypeFile() {
        return myPhenotypeFile.value();
    }

    /**
     * Set Phenotype File. Phenotype File
     *
     * @param value Phenotype File
     *
     * @return this plugin
     */
    public VCAPScanPlugin phenotypeFile(String value) {
        myPhenotypeFile = new PluginParameter<>(myPhenotypeFile, value);
        return this;
    }

    /**
     * Ldak Command
     *
     * @return Ldak Command
     */
    public String ldakCommand() {
        return myLDAKCommand.value();
    }

    /**
     * Set Ldak Command. Ldak Command
     *
     * @param value Ldak Command
     *
     * @return this plugin
     */
    public VCAPScanPlugin ldakCommand(String value) {
        myLDAKCommand = new PluginParameter<>(myLDAKCommand, value);
        return this;
    }

    /**
     * Output Dir
     *
     * @return Output Dir
     */
    public String outputDir() {
        return myOutputDir.value();
    }

    /**
     * Set Output Dir. Output Dir
     *
     * @param value Output Dir
     *
     * @return this plugin
     */
    public VCAPScanPlugin outputDir(String value) {
        myOutputDir = new PluginParameter<>(myOutputDir, value);
        return this;
    }

    /**
     * Whole Matrix
     *
     * @return Whole Matrix
     */
    public String wholeMatrix() {
        return myWholeMatrix.value();
    }

    /**
     * Set Whole Matrix. Whole Matrix
     *
     * @param value Whole Matrix
     *
     * @return this plugin
     */
    public VCAPScanPlugin wholeMatrix(String value) {
        myWholeMatrix = new PluginParameter<>(myWholeMatrix, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "VCAP Scan";
    }

    @Override
    public String getToolTipText() {
        return "VCAP Scan";
    }

}
