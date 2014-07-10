/*
 * ReImputeUpdatedTaxaByFILLINPlugin
 */
package net.maizegenetics.analysis.imputation;

// standard imports for plugins
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.plugindef.AbstractPlugin;
import org.apache.log4j.Logger;
import net.maizegenetics.plugindef.DataSet;
import javax.swing.*;
import java.awt.*;
import net.maizegenetics.plugindef.PluginParameter;

// imports specifically needed for this plugin
import java.util.ArrayList;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import java.util.List;

/**
 * Compares an unfinished HDF5 file containing raw genotypes to a corresponding 
 * unfinished HDF5 file containing FILLIN-imputed genotypes to find new taxa (or 
 * taxa with additional depth) in the raw geno file, then imputes (or reimputes) 
 * these with FILLIN and adds them to (or replaces them in) the imputed geno file.
 * 
 * This is part of the Automated Production Pipeline.
 * 
 * @author jcg233
 */
public class ReImputeUpdatedTaxaByFILLINPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ReImputeUpdatedTaxaByFILLINPlugin.class);

    private PluginParameter<String> rawHDF5GenotypeFile 
        = new PluginParameter.Builder<>("raw", null, String.class)
        .guiName("Raw HDF5 Genotype File")
        .required(true)
        .inFile()
        .description("Input HDF5 (*.h5) file containing raw (unimputed) genotypes")
        .build();
    private PluginParameter<String> imputedHDF5GenotypeFile 
        = new PluginParameter.Builder<>("imp", null, String.class)
        .guiName("Imputed HDF5 Genotype File")
        .required(true)
        .inFile()
        .description("Target HDF5 (*.h5) file containing imputed genotypes to be updated")
        .build();
    private PluginParameter<String> donorDir 
        = new PluginParameter.Builder<>("d",null,String.class)
        .guiName("Donor Dir")
        .inDir()
        .required(true)
        .description("Directory containing donor haplotype files from output of the FILLINFindHaplotypesPlugin. "
                    +"All files with '.gc' in the filename will be read in, only those with matching sites are used")
        .build();

    // TODO: add all possible FILLINImputationPlugin parameters?  It seems that the default parameters were used for maize.
    
    // global variables
    IHDF5Reader rawGenosReader;
    IHDF5Writer impGenosWriter;
    

    public ReImputeUpdatedTaxaByFILLINPlugin() {
        super(null, false);
    }

    public ReImputeUpdatedTaxaByFILLINPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, false);
    }

    @Override
    public String pluginDescription() {
        return 
            "This plugin " +
            "compares an unfinished HDF5 file containing raw genotypes to a corresponding " +
            "unfinished HDF5 file containing FILLIN-imputed genotypes to find new taxa (or " +
            "taxa with additional depth) in the raw geno file, then imputes (or reimputes) " +
            "these with FILLIN and adds them to (or replaces them in) the imputed geno file."
        ;
    }

    @Override
    public DataSet processData(DataSet input) {
        ReImputeUpdatedTaxaByFILLIN();
        fireProgress(100);
        return null;
    }
    
    private void ReImputeUpdatedTaxaByFILLIN() {
        // open raw and target imputed genos (both are unfinished HDF5 genos)
        openInputHDF5GenoFiles();
        
        // compare taxa (exit if no change)
        TaxaList modifiedTaxa = compareRawAndImputedTaxa();
        if (modifiedTaxa.isEmpty()) {
            myLogger.info("No additional or updated taxa were found in the raw genotype input file.");
            return;
        }
        
        // create temporary input HDF5 file (no depth needed) with taxa subset to feed to the FILLINFindHaplotypesPlugin
        String tempInFile = createTempInputFileForFILLIN(modifiedTaxa);
        
        // run FILLINFindHaplotypesPlugin, producing temporary output HDF5 imputed genotypes
        String tempOutFile = runFILLIN(tempInFile);
        
        // replace taxa & genotypes in target HDF5 imputed genotypes file
        replaceTaxaInImputedFile(tempOutFile);
        
        // delete temporary files
        deleteTemporaryFiles(tempInFile, tempOutFile);
        
    }
    
    private void openInputHDF5GenoFiles() {
        myLogger.info("\nOpening input raw genotypes file:\n  "+rawHDF5GenotypeFile()+"\n");
        rawGenosReader=HDF5Factory.openForReading(rawHDF5GenotypeFile());
        myLogger.info("\nOpening target imputed genotypes file:\n  "+imputedHDF5GenotypeFile()+"\n");
        impGenosWriter=HDF5Factory.open(imputedHDF5GenotypeFile());
    }
    
    private TaxaList compareRawAndImputedTaxa() {
        StringBuilder modifiedTaxaReport = new StringBuilder("Modified taxa:\n");
        ArrayList<Taxon> modifiedTaxa = new ArrayList();
        // compare taxa & add to modified taxa if new or changed
        List<String> rawTaxaNames = HDF5Utils.getAllTaxaNames(rawGenosReader);
        for (String taxonName : rawTaxaNames) {
            if (!HDF5Utils.doTaxonCallsExist(impGenosWriter, taxonName)) {
                modifiedTaxa.add(HDF5Utils.getTaxon(rawGenosReader, taxonName));
                modifiedTaxaReport.append("  "+taxonName+" (new taxon)\n");
            } else if (flowcellLaneAdded(taxonName)) {
                modifiedTaxa.add(HDF5Utils.getTaxon(rawGenosReader, taxonName));
                modifiedTaxaReport.append("  "+taxonName+" (additional depth)\n");
            }
        }
        if (!modifiedTaxa.isEmpty()) myLogger.info(modifiedTaxaReport.toString());
        return new TaxaListBuilder().addAll(modifiedTaxa).sortTaxaAlphabetically().build();
    }
    
    private boolean flowcellLaneAdded(String taxonName) {
        Taxon rawTaxon = HDF5Utils.getTaxon(rawGenosReader, taxonName);
        if (rawTaxon == null) {
            throw new IllegalStateException("No corresponding Taxon found in the raw genotype file for the existing taxon name: "+taxonName);
        }
        Taxon impTaxon = HDF5Utils.getTaxon(impGenosWriter, taxonName);
        if (impTaxon == null) return true;
        String[] rawFlowCellLanes = (String[]) rawTaxon.getAnnotation("Flowcell_Lane");
        String[] impFlowCellLanes = (String[]) impTaxon.getAnnotation("Flowcell_Lane");
        for (String rawFlowCellLane : rawFlowCellLanes) {
            boolean found = false;
            for (String impFlowCellLane : impFlowCellLanes) {
                if(impFlowCellLane.equals(rawFlowCellLane)) {
                    found = true;
                    continue;
                }
            }
            if (!found) return true;
        }
        return false;
    }
    
    private String createTempInputFileForFILLIN(TaxaList modifiedTaxa) {
        return null;
    }
    
    private String runFILLIN(String tempInFile) {
        return null;
    }
    
    private void replaceTaxaInImputedFile(String tempOutFile) {
        ;
    }
    
    private void deleteTemporaryFiles(String tempInFile, String tempOutFile) {
        ;
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Update imputed genotypes";
    }

    @Override
    public String getToolTipText() {
        return "Update imputed genotypes file based on modified/new taxa in raw genotypes file";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ReImputeUpdatedTaxaByFILLINPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public DataSet runPlugin(DataSet input) {
        return (DataSet) performFunction(input).getData(0).getData();
    }

    /**
     * Input HDF5 (*.h5) file containing raw (unimputed) genotypes
     *
     * @return Raw HDF5 Genotype File
     */
    public String rawHDF5GenotypeFile() {
        return rawHDF5GenotypeFile.value();
    }

    /**
     * Set Raw HDF5 Genotype File. Input HDF5 (*.h5) file
     * containing raw (unimputed) genotypes
     *
     * @param value Raw HDF5 Genotype File
     *
     * @return this plugin
     */
    public ReImputeUpdatedTaxaByFILLINPlugin rawHDF5GenotypeFile(String value) {
        rawHDF5GenotypeFile = new PluginParameter<>(rawHDF5GenotypeFile, value);
        return this;
    }

    /**
     * Target HDF5 (*.h5) file containing imputed genotypes
     * to be updated
     *
     * @return Imputed HDF5 Genotype File
     */
    public String imputedHDF5GenotypeFile() {
        return imputedHDF5GenotypeFile.value();
    }

    /**
     * Set Imputed HDF5 Genotype File. Target HDF5 (*.h5)
     * file containing imputed genotypes to be updated
     *
     * @param value Imputed HDF5 Genotype File
     *
     * @return this plugin
     */
    public ReImputeUpdatedTaxaByFILLINPlugin imputedHDF5GenotypeFile(String value) {
        imputedHDF5GenotypeFile = new PluginParameter<>(imputedHDF5GenotypeFile, value);
        return this;
    }

    /**
     * Directory containing donor haplotype files from output
     * of the FILLINFindHaplotypesPlugin. All files with '.gc'
     * in the filename will be read in, only those with matching
     * sites are used
     *
     * @return Donor Dir
     */
    public String donorDir() {
        return donorDir.value();
    }

    /**
     * Set Donor Dir. Directory containing donor haplotype
     * files from output of the FILLINFindHaplotypesPlugin.
     * All files with '.gc' in the filename will be read in,
     * only those with matching sites are used
     *
     * @param value Donor Dir
     *
     * @return this plugin
     */
    public ReImputeUpdatedTaxaByFILLINPlugin donorDir(String value) {
        donorDir = new PluginParameter<>(donorDir, value);
        return this;
    }

}
