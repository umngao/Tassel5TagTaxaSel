package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalTransformPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.LoggingUtils;

public class TransformDataPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(TransformDataPlugin.class);

	public DataSet processData(DataSet input){
		if (input.getSize() > 1) {
			throw new IllegalArgumentException("Select a single dataset for transformation.");
		}
		
		List<Datum> myData = input.getDataOfType(GenotypeTable.class);
		
		if (myData.size() == 1) {
			GenotypeTable myGenotype = (GenotypeTable) myData.get(0).getData();
			//TODO call GenotypeTransformPlugin
			
		}
		
		myData = input.getDataOfType(Phenotype.class);
		if (myData.size() == 1) {
			Phenotype myPhenotype = (Phenotype) myData.get(0).getData();
			//TODO call PhenotypeTransformPlugin
			
		}
		
		return null;
	}
	
	public TransformDataPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}

	@Override
	public ImageIcon getIcon() {
        URL imageURL = NumericalTransformPlugin.class.getResource("Transform.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL); 
        }
	}

	@Override
	public String getButtonName() {
		return "Transform";
	}

	@Override
	public String getToolTipText() {
		return "Transform phenotypes or convert genotypes to probabilities";
	}

}
