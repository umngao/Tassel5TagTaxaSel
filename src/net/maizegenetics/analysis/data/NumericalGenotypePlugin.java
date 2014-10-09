/*
* Plugin for numerical genotype.
* @author - Janu Verma
* jv367@cornell.edu
 */

package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

import java.awt.Frame;
import java.lang.Override;
import java.lang.String;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;



public class NumericalGenotypePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(NumericalGenotypePlugin.class);

    /** Creates a new instance of the NumericalGenotypePlugin */
    public NumericalGenotypePlugin(Frame parentFrame, boolean isInteractive){
        super(parentFrame, isInteractive);
    }


    @Override
    public DataSet processData(DataSet input) {

        List<Datum> datumList = input.getDataOfType(GenotypeTable.class);

        //check size of datumList, throw error if not equal to one
        if (datumList.size() != 1){
            throw new IllegalArgumentException("NumericalGenotypePlugin: select exactly one input dataset.");
        }

        //load the GenotypeTable.
        GenotypeTable myGenotype = (GenotypeTable) datumList.get(0).getData();

        int nsites = myGenotype.numberOfSites();
        int ntaxa = myGenotype.numberOfTaxa();

        ReferenceProbability myProb;
        double[][] data = null;
        //if GenotypeTable has no ReferenceProbablity, compute the probabilities.
        if (!myGenotype.hasReferenceProbablity()) {
            data = GenotypeTableUtils.convertGenotypeToDoubleProbability(myGenotype, true);
        }

        //build new ReferenceProbability
        ReferenceProbabilityBuilder refBuilder = ReferenceProbabilityBuilder.getInstance(ntaxa, nsites, myGenotype.taxa());
        for (int t = 0; t < ntaxa; t++) {
            float[] values = new float[nsites];

            for (int s = 0; s < nsites; s++) {
                values[s] = (float) data[s][t];
            }
            refBuilder.addTaxon(t, values);
        }

        //build new GenotypeTable
        GenotypeTable numGenotype = GenotypeTableBuilder.getInstance(myGenotype.genotypeMatrix(), myGenotype.positions(),
                myGenotype.taxa(), myGenotype.depth(), myGenotype.alleleProbability(), refBuilder.build(), myGenotype.dosage(),
                myGenotype.annotations());

        String name = datumList.get(0).getName() + "_with_Probability";
        String comment = "A Reference Probability has been computed.";
        Datum newDatum = new Datum(name, numGenotype, comment);
        return new DataSet(newDatum, this);

    }


    @Override
    public String pluginDescription(){
        return "This plugin creates a  numerical genotype table for input genotype data";
    }

    @Override
    public String getCitation(){
        return null;
    }

    @Override
    public  ImageIcon getIcon(){
        return null;
    }

    @Override
    public String getButtonName(){
        return "Numerical Genotype";
    }


    @Override
    public String getToolTipText(){
        return "Numerical Genotype";
    }

}



