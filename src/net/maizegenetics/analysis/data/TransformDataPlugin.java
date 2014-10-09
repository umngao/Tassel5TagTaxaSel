package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.List;

import javax.swing.ImageIcon;


import net.maizegenetics.analysis.numericaltransform.NumericalTransformPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class TransformDataPlugin extends AbstractPlugin {

	public TransformDataPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}

	public DataSet processData(DataSet input){
		if (input.getSize() > 1) {
			throw new IllegalArgumentException("TransformDataPlugin: Select a single dataset for transformation.");
		}

		List<Datum> myData = input.getDataOfType(GenotypeTable.class);
		if (myData.size() == 1) {
			NumericalGenotypePlugin ngp = new NumericalGenotypePlugin(getParentFrame(), isInteractive());
			return ngp.processData(input);
		}

		myData = input.getDataOfType(Phenotype.class);
		if (myData.size() == 1) {
			PhenotypeTransformPlugin ptp = new PhenotypeTransformPlugin(getParentFrame(), isInteractive());
			return ptp.processData(input);
		}

		throw new IllegalArgumentException("TransformDataPlugin: the dataset selected is of the wrong type.");
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
