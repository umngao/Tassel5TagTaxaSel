/**
 * 
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.util.ArrayList;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.log4j.Logger;

/**
 * Plugin to thin out sites based on their physical position on the chromosome.
 * Expects the user to enter a minimum distance between sites in base pair.
 * Accepts and writes .hmp.txt, .hmp.txt.gz, .vcf, .vcf.gz and .h5 files.
 * 
 * @author StefanReuscher
 * @author JeffGlaubitz
 * @author lcj34
 *
 */
public class ThinSitesByPositionPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(GenosToABHPlugin.class);
    private PluginParameter<String> infile = new PluginParameter.Builder<>("i", null, String.class)
            .required(true)
            .inFile()
            .guiName("Input file")
            .description("Input genotype file to be thinned").build();
    private PluginParameter<String> outfile= new PluginParameter.Builder<>("o", null, String.class)
            .required(true)
            .outFile()
            .guiName("Output file")
            .description("Output genotype file").build();
    private PluginParameter<Integer> minDist= new PluginParameter.Builder<>("minDist", null, Integer.class)
            .required(true)
            .guiName("Minimum distance")
            .description("Minimum distance in bp between adjacent sites")
            .build();
    
    public ThinSitesByPositionPlugin() {
        super(null, false);
    }

    public ThinSitesByPositionPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet processData(DataSet input) {
        try {
            //Read the genotypes file - create genotype table
            GenotypeTable genos= ImportUtils.readGuessFormat(infile());
            String[] sitesToKeep = findSiteNamesToKeep(genos);         

            GenotypeTable filteredGenos = null;
            if ((sitesToKeep != null) && (sitesToKeep.length != 0)) {
                filteredGenos = FilterGenotypeTable.getInstance(genos, sitesToKeep);
                writeOutputGenos(filteredGenos);
            } else {
                myLogger.warn("WARNING - no sites kept, no output file written !!");
            }
                                         
        } finally {
            fireProgress(100);
        }       
        return null;
    }
    

    private String[] findSiteNamesToKeep(GenotypeTable genos){
        // This method traverses the Positions in each chromosome,
        // adding to a list those positions within the chromosome
        // that are the user specified distance apart.
        ArrayList<String> sitesToKeep = new ArrayList<String>();
        PositionList positions = genos.positions();
        int chrom = -1;
        int prevPos = -1;
        for (Position currPos : positions){
            if (currPos.getChromosome().getChromosomeNumber() != chrom) {
                // Always keep the first position in a chromosome
                sitesToKeep.add(currPos.getSNPID());
                chrom = currPos.getChromosome().getChromosomeNumber();
                prevPos = currPos.getPosition();
            } else {
                if (currPos.getPosition() - prevPos >= minDist()) {;
                    sitesToKeep.add(currPos.getSNPID());
                    prevPos = currPos.getPosition();
                } // do nothing if distance between positions < minDistß
            }           
        }  
        myLogger.info("ThinSitesByPosition: original number of sites: " + positions.size() 
                + " number of sites kept: " + sitesToKeep.size());
        String[] stringArray = sitesToKeep.toArray(new String[sitesToKeep.size()]);
        return stringArray;
    }
    
    private void writeOutputGenos (GenotypeTable genos){
        // Output format is based on name given for output file
        // If the extension is not a TASSEL supported extension, default to hapmap
        String fileName = outFile();
        if (fileName.endsWith(".h5")) {
            ExportUtils.writeGenotypeHDF5(genos,fileName,true);            
        } else if (fileName.endsWith("hmp.txt.gz") || fileName.endsWith("hmp.txt")) {
            ExportUtils.writeToHapmap(genos, fileName);
        } else if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
            ExportUtils.writeToVCF(genos, fileName, true);
        } else {
            fileName = fileName + ".hmp.txt";
            myLogger.warn("File extension not recognized, writing to hapmap file:\n    " + fileName);
            ExportUtils.writeToHapmap(genos, fileName);
        }
    }
    
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ThinSitesByPositionPlugin.class);
    // }


    /**
     * Input genotype file to be thinned
     *
     * @return Input file
     */
    public String infile() {
        return infile.value();
    }

    /**
     * Set Input file. Input genotype file to be thinned
     *
     * @param value Input file
     *
     * @return this plugin
     */
    public ThinSitesByPositionPlugin infile(String value) {
        infile = new PluginParameter<>(infile, value);
        return this;
    }

    /**
     * Output genotype file
     *
     * @return Output file
     */
    public String outFile() {
        return outfile.value();
    }

    /**
     * Set Output file. Output genotype file
     *
     * @param value Output file
     *
     * @return this plugin
     */
    public ThinSitesByPositionPlugin outfile(String value) {
        outfile = new PluginParameter<>(outfile, value);
        return this;
    }

    /**
     * Minimum distance in bp between adjacent sites
     *
     * @return Minimum distance
     */
    public Integer minDist() {
        return minDist.value();
    }

    /**
     * Set Minimum distance. Minimum distance in bp between
     * adjacent sites
     *
     * @param value Minimum distance
     *
     * @return this plugin
     */
    public ThinSitesByPositionPlugin minDist(Integer value) {
        minDist = new PluginParameter<>(minDist, value);
        return this;
    }

    @Override
    public String getCitation() {
        return "Stefan Reuscher, Jeff Glaubitz, Lynn Johnson (2015) First Annual TASSEL Hackathon";
    }

//     The following getters and setters were auto-generated.
//     Please use this method to re-generate.
    
//     public static void main(String[] args) {
//         GeneratePluginCode.generate(ThinSitesByPositionPlugin.class);
//     }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
//    public void runPlugin(DataSet input) {
//        return (void) performFunction(input).getData(0).getData();
//    }
}
