package net.maizegenetics.analysis.data;

// standard imports for plugins
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;
import javax.swing.*;
import java.awt.*;
import net.maizegenetics.plugindef.PluginParameter;
//import net.maizegenetics.plugindef.GeneratePluginCode;

// specifically needed for this plugin
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;

/**
 * Opens an "unfinished" HDF5 genos file, on which .closeUnfinished() was called rather than .build(),
 * makes a copy of it and then finalizes it by calling .build().
 * 
 * If provided, root level annotations DataSetName and DataSetDescription are also added.
 *
 * @author Jeff Glaubitz (jcg233@cornell.edu)
 */
public class BuildUnfinishedHDF5GenotypesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(BuildUnfinishedHDF5GenotypesPlugin.class);

    private PluginParameter<String> inputGenotypes = new PluginParameter.Builder<>("i", null, String.class)
        .guiName("Input file")
        .required(true)
        .inFile()
        .description("Input, unfinished HDF5 genotype (*.h5) file to be fininalized")
        .build();
    private PluginParameter<String> outputGenotypes = new PluginParameter.Builder<>("o", null, String.class)
        .guiName("Output file")
        .required(true)
        .outFile()
        .description("Output, finished HDF5 genotype (*.h5) file which can be opened with the TASSEL5 GUI")
        .build();
    private PluginParameter<String> dataSetName = new PluginParameter.Builder<>("name", null, String.class)
        .guiName("Data set name")
        .required(false)
        .description("(Optional) Short data set name to be added as an root level annotation under \"/DataSetName\"")
        .build();
    private PluginParameter<String> dataSetDescription = new PluginParameter.Builder<>("desc", null, String.class)
        .guiName("Data set description")
        .required(false)
        .description("(Optional) Short data set description to be added as an root level annotation under \"/DataSetDescription\"")
        .build();

    public BuildUnfinishedHDF5GenotypesPlugin() {
        super(null, false);
    }

    public BuildUnfinishedHDF5GenotypesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, false);
    }

    @Override
    public DataSet processData(DataSet input) {
        String message = buildUnfinishedHDF5Genotypes();
        if(message != null) {
            myLogger.error(message);
            try {Thread.sleep(500);} catch(Exception e) {}
            throw new IllegalStateException(message);
        }
        fireProgress(100);
        return null;
    }

    private String buildUnfinishedHDF5Genotypes() {
        myLogger.info("\n\nBuildUnfinishedHDF5GenotypesPlugin:\nFinalizing the following HDF5 genotype file:\n  "
                +inputFile()+"\n\n");
        IHDF5Reader h5Reader = HDF5Factory.openForReading(inputFile());
        GenotypeTableBuilder genoTable = GenotypeTableBuilder.getTaxaIncremental(PositionListBuilder.getInstance(h5Reader), outputFile());
        TaxaList taxa = new TaxaListBuilder().buildFromHDF5Genotypes(h5Reader);
        for (Taxon taxon : taxa) {
            byte[] genos = HDF5Utils.getHDF5GenotypesCalls(h5Reader, taxon.getName());
            byte[][] depth = HDF5Utils.getHDF5GenotypesDepth(h5Reader, taxon.getName());
            genoTable.addTaxon(taxon, genos, depth);
        }
        
        // need to add these via GenotypeTableBuilder
        if (dataSetName() != null) {
//            HDF5Utils.writeHDF5DataSetName(null, null);
        }
        if (dataSetDescription() != null) {
            ;
        }
        
        genoTable.build();
        myLogger.info("\n\nBuildUnfinishedHDF5GenotypesPlugin: Finished finalizing HDF5 genotypes in the following output file:\n  "
                +outputFile()+"\n\n");
        return null;
    }
    
    private String parseDataSetName(String dataSetName) {
        return null;
    }
    
    private String parseDataSetDescription(String dataSetDescrip, int nSNPs, int nTaxa) {
        return null;
    }

    @Override
    public ImageIcon getIcon() {
//        URL imageURL = SeparatePlugin.class.getResource("/net/maizegenetics/analysis/images/Merge.gif");
//        if (imageURL == null) {
            return null;
//        } else {
//            return new ImageIcon(imageURL);
//        }
    }

    @Override
    public String getButtonName() {
        return "Split chromosomes from HDF5 genotype file";
    }

    @Override
    public String getToolTipText() {
        return "Split chromosomes from HDF5 genotype file";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(BuildUnfinishedHDF5GenotypesPlugin.class);
    // }

    /**
     * Input, unfinished HDF5 genotype (*.h5) file to be fininalized
     *
     * @return Input file
     */
    public String inputFile() {
        return inputGenotypes.value();
    }

    /**
     * Set Input file. Input, unfinished HDF5 genotype (*.h5)
     * file to be fininalized
     *
     * @param value Input file
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin inputFile(String value) {
        inputGenotypes = new PluginParameter<>(inputGenotypes, value);
        return this;
    }

    /**
     * Output, finished HDF5 genotype (*.h5) file which can
     * be opened with the TASSEL5 GUI
     *
     * @return Output file
     */
    public String outputFile() {
        return outputGenotypes.value();
    }

    /**
     * Set Output file. Output, finished HDF5 genotype (*.h5)
     * file which can be opened with the TASSEL5 GUI
     *
     * @param value Output file
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin outputFile(String value) {
        outputGenotypes = new PluginParameter<>(outputGenotypes, value);
        return this;
    }

    /**
     * (Optional) Short data set name to be added as an root
     * level annotation under "/DataSetName"
     *
     * @return Data set name
     */
    public String dataSetName() {
        return dataSetName.value();
    }

    /**
     * Set Data set name. (Optional) Short data set name to
     * be added as an root level annotation under "/DataSetName"
     *
     * @param value Data set name
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin dataSetName(String value) {
        dataSetName = new PluginParameter<>(dataSetName, value);
        return this;
    }

    /**
     * (Optional) Short data set description to be added as
     * an root level annotation under "/DataSetDescription"
     *
     * @return Data set description
     */
    public String dataSetDescription() {
        return dataSetDescription.value();
    }

    /**
     * Set Data set description. (Optional) Short data set
     * description to be added as an root level annotation
     * under "/DataSetDescription"
     *
     * @param value Data set description
     *
     * @return this plugin
     */
    public BuildUnfinishedHDF5GenotypesPlugin dataSetDescription(String value) {
        dataSetDescription = new PluginParameter<>(dataSetDescription, value);
        return this;
    }

}
