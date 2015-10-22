/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.CoreGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.MaskGenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;

/**
 *
 * @author jcg233
 */
public class SetLowDepthGenosToMissingPlugin extends net.maizegenetics.plugindef.AbstractPlugin {

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
    private String inputGenosName = null;
    
    public SetLowDepthGenosToMissingPlugin(Frame parentFrame, boolean isInteractive) { 
        super(parentFrame, isInteractive);
    } 

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException(
                "SetLowDepathGenotypesToMissingPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);
        
        inputGenosName = genotypeTables.get(0).getName();

        System.out.println("Input genotype name: " + inputGenosName);
        
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
    
    @Override
    public DataSet processData(DataSet input) {
        int numberOfSites = inputGenotypes.numberOfSites();
        MaskGenotypeTableBuilder mgtb = new MaskGenotypeTableBuilder(inputGenotypes);


        for (int siteIndex = 0; siteIndex < numberOfSites; siteIndex++) {
            int[][] allelesByFreq = inputGenotypes.allelesSortedByFrequency(siteIndex);
            int numAlleles = allelesByFreq[0].length;
            int[] allelesAtSite = new int[numAlleles];
            
            for (int alleleIndex = 0; alleleIndex < numAlleles; alleleIndex++){
                allelesAtSite[alleleIndex] = allelesByFreq[0][alleleIndex];
            }
            
            for (Taxon inTaxon : inputGenotypes.taxa()) {
                int taxonIndex = inputGenotypes.taxa().indexOf(inTaxon);
                int[] alleleDepths = inputGenotypes.depthForAlleles(taxonIndex, siteIndex);
                int depthSum = 0;
                for ( int allele : allelesAtSite){
                   depthSum += alleleDepths[allele];
                }
                if (depthSum < minDepth()) {
                    mgtb.set(taxonIndex, siteIndex);
                }
            }
        }

        String outGenosName = inputGenosName + "MinDepth" + minDepth();
        GenotypeTable outGenos = mgtb.build();
        return new DataSet(new Datum(outGenosName, outGenos, null), null);
    }
    
    @Override
    public ImageIcon getIcon() {
        URL imageURL = SetLowDepthGenosToMissingPlugin.class
            .getResource("/net/maizegenetics/analysis/images/lowDepthToMissing.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "SetLowDepthGenosToMissing";
    }

    @Override
    public String getToolTipText() {
        return "Set genotypes below a minimum depth to missing";
    }
 
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(SetLowDepthGenosToMissingPlugin.class);
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
    public SetLowDepthGenosToMissingPlugin minDepth(Integer value) {
        minDepth = new PluginParameter<>(minDepth, value);
        return this;
    }
    
}
