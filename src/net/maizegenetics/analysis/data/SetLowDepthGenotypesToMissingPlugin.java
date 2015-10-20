/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

/**
 *
 * @author jcg233
 */
public class SetLowDepthGenotypesToMissingPlugin extends net.maizegenetics.plugindef.AbstractPlugin {

    @Override
    public String pluginDescription() {
        return "Sets each genotype in the input genotypes to missing if "
            +  "the underlying allelic depth is below a specified minimum.";
    }
    
    private PluginParameter<Integer> minDepth = 
        new PluginParameter.Builder<>("minDepth",null,Integer.class)
            .guiName("minimum genotype depth")
            .required(true)
            .description("Minimum depth, below which genotypes are set to missing")
            .build();
    
    private GenotypeTable inputGenotypes = null;
    
    public SetLowDepthGenotypesToMissingPlugin(Frame parentFrame, boolean isInteractive) { 
        super(parentFrame, isInteractive);
    } 

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException(
                "SetLowDepathGenotypesToMissingPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);
        if (genotypeTables.size() == 1) {
            inputGenotypes
                = (GenotypeTable) genotypeTables.get(0).getData();
            if (! inputGenotypes.hasDepth()){
                throw new IllegalArgumentException("SetLowDepthGenotypesToMissingPlugin: preProcessParameters: Please select a Genotype Table with allele depth information.");
            }
        } else {
            throw new IllegalArgumentException(
                "SetLowDepthGenotypesToMissingPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }
    
    public DataSet processData(DataSet input) {
        for (int taxonIndex = 0; taxonIndex < inputGenotypes.numberOfTaxa(); taxonIndex++) {
            for (int siteIndex = 0; siteIndex < inputGenotypes.numberOfSites(); siteIndex++) {
                int[][] allelesByFreq = inputGenotypes.allelesSortedByFrequency(siteIndex);
                int[] alleleDepths = inputGenotypes.depthForAlleles(taxonIndex, siteIndex);
                int depthSum = 0;
                for ( int allele = 0; allele < alleleDepths.length; allele++){
                    depthSum += alleleDepths[allele];
                }
                System.out.println("Allele depth sum: " + depthSum);
            }
            
        }
        return null;
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "SetLowDepthGenosToMissing";
    }

    @Override
    public String getToolTipText() {
        return "Set Genotypes below a minimum depth to missing";
    }
 
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(SetLowDepthGenotypesToMissingPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Min Depth
     *
     * @return Min Depth
     */
    public Integer minDepth() {
        return minDepth.value();
    }

    /**
     * Set Min Depth. Min Depth
     *
     * @param value Min Depth
     *
     * @return this plugin
     */
    public SetLowDepthGenotypesToMissingPlugin minDepth(Integer value) {
        minDepth = new PluginParameter<>(minDepth, value);
        return this;
    }
    
}
