package net.maizegenetics.analysis.numericaltransform;

import java.awt.Frame;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

public class AvgPhenotypeByTaxaPlugin extends AbstractPlugin{
    
    private PluginParameter<Boolean> addSmallValue =
            new PluginParameter.Builder<>("addSmallValue", false, Boolean.class)
            .build();
    
    public AvgPhenotypeByTaxaPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    public DataSet processData(DataSet input) {
        //do whatever your plugin does
        List<Datum> datumList = input.getDataOfType(Phenotype.class);

        //check size of datumList, throw error if not equal to one
        if (datumList.size() != 1){
            throw new IllegalArgumentException("AvgPhenotypeByTaxaPlugin: select exactly one phenotype dataset to average.");
        }
        
        AvgPhenotype avgPheno = new AvgPhenotype();
        try {
        Phenotype myPhenotype = avgPheno.averagePheno((Phenotype) datumList.get(0).getData(),addSmallValue.value());
        
        if (myPhenotype != null) {
            String name = "PhenotypeTransformed";

            Datum td = new Datum(name, myPhenotype, null);
            
            //todo need to add logic of directories.
            DataSet tds = new DataSet(td, this);
            return tds;
        }
        }
        catch(Exception e) {
            throw new IllegalArgumentException("Error building Phenotype Object.  "
                    + "Check to make sure the dimensions of your data match the number of taxa and the number of attributes.");
        }
        return null;  // Note: this can return null
        //return myPhenotype;  // Note: this can return null
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }
    @Override
    public String getButtonName() {
        return "Average Phenotypes By Taxa";
    }
    @Override
    public String getToolTipText() {
        return "Averages Phenotypes for Duplicate Taxa.";
    }
    
    @Override
    public String getCitation() {
        return "Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. "
                + "(2007) TASSEL: Software for association mapping of complex traits in diverse "
                + "samples. Bioinformatics 23:2633­2635.";
    }
    
    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel­5­source/wiki/UserManual/Kinship/Kinship"; 
    } 
    
 // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(AvgPhenotypeByTaxaPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public Phenotype runPlugin(DataSet input) {
        return (Phenotype) performFunction(input).getData(0).getData();
    }

    /**
     * Add Small Value
     *
     * @return Add Small Value
     */
    public Boolean addSmallValue() {
        return addSmallValue.value();
    }

    /**
     * Set Add Small Value. Add Small Value
     *
     * @param value Add Small Value
     *
     * @return this plugin
     */
    public AvgPhenotypeByTaxaPlugin addSmallValue(Boolean value) {
        addSmallValue = new PluginParameter<>(addSmallValue, value);
        return this;
    }
}
