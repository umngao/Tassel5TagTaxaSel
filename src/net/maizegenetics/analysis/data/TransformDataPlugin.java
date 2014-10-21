package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.stream.Stream;

import javax.swing.ImageIcon;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.analysis.numericaltransform.NumericalTransformPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.OpenBitSet;

public class TransformDataPlugin extends AbstractPlugin {
	private enum BASE {natural, base_2, base_10};
	
	private double power = 1;
	private List<CategoricalAttribute> byFactor;
	private BASE myBase = BASE.natural;
	
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
//			PhenotypeTransformPlugin ptp = new PhenotypeTransformPlugin(getParentFrame(), isInteractive());
//			return ptp.processData(input);
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

	//method implementing transformation using stream method
	public static NumericAttribute transformUsingStream(NumericAttribute original, DoubleUnaryOperator transformOp) {
		double[] values = original.stream().map(transformOp).toArray();
		return new NumericAttribute(original.name(), AssociationUtils.convertDoubleArrayToFloat(values), original.missing());
	}
	
	public NumericAttribute powerTransform(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] transValues = new float[n];
		
		for (int i = 0; i < n; i++) transValues[i] = (float) Math.pow(originalValues[i], power);
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public NumericAttribute logTransform(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] transValues = new float[n];
		double divisor;
		
		switch (myBase) {
		case base_10:
			divisor = Math.log(10);
			break;
		case base_2:
			divisor = Math.log(2);
			break;
		default:
			divisor = 1;
		}
		
		for (int i = 0; i < n; i++) {
			switch (myBase) {
			case natural:
				transValues[i] = (float) Math.log(originalValues[i]);
				break;
			case base_10:
			case base_2:
				transValues[i] = (float) (Math.log(originalValues[i]) / divisor);
				break;
			}
		}
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public NumericAttribute lnTransform(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] transValues = new float[n];
		
		for (int i = 0; i < n; i++) transValues[i] = (float) Math.log(originalValues[i]);
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}

	public NumericAttribute standardize(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] meanSD = meanStdDev(originalValues);
		float[] transValues = new float[n];
		
		for (int i = 0; i < n; i++) transValues[i] = (transValues[i] - meanSD[0])/meanSD[1];
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public float[] meanStdDev(float[] data) {
		int n = data.length;
		float sum = 0;
		float sumsq = 0;
		int notMissingCount = 0;
		for (int i = 0; i < n; i++) {
			float val = data[i];
			if (!Float.isNaN(val)) {
				sum += val;
				sumsq += val * val;
				notMissingCount++;
			}
			
		}
		float mean = sum / notMissingCount;
		float sdev = (float) Math.sqrt(sumsq - sum / notMissingCount * sum);
		return new float[]{mean, sdev};
	}
	
	public NumericAttribute standardize(NumericAttribute original, List<PhenotypeAttribute> byFactors) {
		
		List<int[]> subsetList = subsets(byFactors);
		float[] stdData = Arrays.copyOf(original.floatValues(), original.size());
		for (int[] subset : subsetList) {
			int n = subset.length;
			float[] subsetData = new float[n];
			for (int i = 0; i < n; i++) subsetData[i] = stdData[subset[i]];
			float[] meanSD = meanStdDev(stdData);
			for (int i = 0; i < n; i++) stdData[subset[i]] = (stdData[subset[i]] - meanSD[0]) / meanSD[1] ;
		}
		
		return new NumericAttribute(original.name(), stdData, original.missing());
	}
	
	public List<int[]> subsets(List<PhenotypeAttribute> byFactors) {
		class subset {
			int[] levels;
			subset(int[] levels) { this.levels = levels; }
			public boolean equals(Object other) {
				if (other instanceof int[]) return Arrays.equals(levels, (int[]) other);
				return false;
			}
			
			@Override
			public int hashCode() {
				int hc = 0;
				int mult = 1;
				for (int i : levels) {
					hc += mult * Integer.hashCode(i);
					mult *= 10;
				}
				return hc;
			}
		}
		
		ArrayList<CategoricalAttribute> caList = new ArrayList<>();
		for (PhenotypeAttribute pa : byFactors) caList.add((CategoricalAttribute) pa);
		int nobs = byFactors.get(0).size();
		OpenBitSet missing = new OpenBitSet(nobs);
		for (PhenotypeAttribute pa : byFactors) missing.or(pa.missing());
		
		int nfactors = byFactors.size();
		Multimap<subset, Integer> subsetMap = HashMultimap.create();
		for (int obs = 0; obs < nobs; obs++) if (!missing.fastGet(obs)) {
			int[] levels = new int[nfactors];
			int count = 0;
			for (CategoricalAttribute ca : caList) {
				levels[count++] = ca.intValue(obs);
			}
			subsetMap.put(new subset(levels), obs);
		}
		
		ArrayList<int[]> subsetList = new ArrayList<>();
		for (subset sub : subsetMap.keySet()) {
			subsetList.add(subsetMap.get(sub)
					.stream()
					.mapToInt(Integer::intValue)
					.toArray()); 
		}
		return subsetList;
	}

}
