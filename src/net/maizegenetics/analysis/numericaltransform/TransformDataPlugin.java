package net.maizegenetics.analysis.numericaltransform;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.swing.ImageIcon;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.OpenBitSet;

public class TransformDataPlugin extends AbstractPlugin {
	public enum BASE {natural, base_2, base_10};
	
	private List<NumericAttribute> traitsToTransform;
	private List<CategoricalAttribute> byFactor;
	private boolean logTransform = false;
	private boolean powerTransform = false;
	private boolean standardize = false;
	private BASE myBase = BASE.natural;
	private double power = 1;
	private static final double log2 = Math.log(2);
		
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
			Phenotype myPhenotype = (Phenotype) myData.get(0).getData();
			if (isInteractive()) {
				
				List<NumericAttribute> numericAttributes = Stream.concat(myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).stream(), 
						myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor).stream())
						.map(pa -> (NumericAttribute) pa)
						.collect(Collectors.toList());
				
				List<CategoricalAttribute> catAttributes = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor).stream()
						.map(pa -> (CategoricalAttribute) pa)
						.collect(Collectors.toList());
				
				TransformDataDialog tdd = new TransformDataDialog(getParentFrame(), numericAttributes, catAttributes);
				tdd.setVisible(true);

				traitsToTransform = tdd.traitsToTransform();
				byFactor = tdd.factorsForStandardizing();
				logTransform = tdd.logTransformation();
				powerTransform = tdd.powerTransformation();
				standardize = tdd.standardize();
				myBase = tdd.base();
				power = tdd.exponent();
			}
			if (logTransform || powerTransform || standardize) return transformTraits(myPhenotype, myData.get(0));
			else return null;
			
		}

		throw new IllegalArgumentException("TransformDataPlugin: the dataset selected is of the wrong type.");
	}

	public DataSet transformTraits(Phenotype myPhenotype, Datum myData) {
		//use a sequential stream, because the order of the attributes needs to stay the same
		List<PhenotypeAttribute> myNewAttributes = myPhenotype.attributeListCopy().stream()
				.map(a -> transformAttribute(a))
				.collect(Collectors.toList());
		
		Phenotype transformedPhenotype = new PhenotypeBuilder().fromAttributeList(myNewAttributes, myPhenotype.typeListCopy()).build().get(0);
		
		StringBuilder nameBuilder = new StringBuilder();
		nameBuilder.append("transformed_").append(myData.getName());
		
		StringBuilder commentBuilder = new StringBuilder();
		commentBuilder.append("Phenotypes transformed from ");
		commentBuilder.append(myData.getName()).append("\n");
//		commentBuilder.append(myData.getComment());
		commentBuilder.append("The following traits were transformed by ");
		if (powerTransform) commentBuilder.append("using a power ").append(power).append(" transformation:\n");
		else if (logTransform) commentBuilder.append("using a ").append(myBase.name()).append(" log transformation:\n");
		if (standardize) commentBuilder.append("standardizing.\n");
		for (NumericAttribute na : traitsToTransform) commentBuilder.append(na.name()).append("\n");
		
		return new DataSet(new Datum(nameBuilder.toString(), transformedPhenotype, commentBuilder.toString()), this);
	}
	
	public PhenotypeAttribute transformAttribute(PhenotypeAttribute myAttribute) {
		if (!(myAttribute instanceof NumericAttribute)) return myAttribute;
		NumericAttribute myNumericAttribute = (NumericAttribute) myAttribute;
		if (!traitsToTransform.contains(myNumericAttribute)) return myAttribute;
		
		if (powerTransform) myNumericAttribute = powerTransform(myNumericAttribute);
		else if (logTransform) myNumericAttribute = logTransform(myNumericAttribute);
		
		if (standardize) {
			if (byFactor.size() > 0) return standardize(myNumericAttribute, byFactor);
			return standardize(myNumericAttribute);
		}
		
		return myNumericAttribute;
	}
	
	@Override
	public ImageIcon getIcon() {
        URL imageURL = TransformDataPlugin.class.getResource("/net/maizegenetics/analysis/images/Transform.gif");
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
	
	public NumericAttribute standardize(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] meanSD = meanStdDev(originalValues);
		float[] transValues = new float[n];
		
		for (int i = 0; i < n; i++) transValues[i] = (originalValues[i] - meanSD[0])/meanSD[1];
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public float[] meanStdDev(float[] data) {
		int n = data.length;
		double sum = 0;
		double sumsq = 0;
		int notMissingCount = 0;
		for (int i = 0; i < n; i++) {
			double val = data[i];
			if (!Double.isNaN(val)) {
				sum += val;
				sumsq += val * val;
				notMissingCount++;
			}
			
		}
		float mean = (float) sum / notMissingCount;
		float sdev = (float) Math.sqrt((sumsq - sum / notMissingCount * sum) / (notMissingCount - 1));
		return new float[]{mean, sdev};
	}
	
	public NumericAttribute standardize(NumericAttribute original, List<CategoricalAttribute> byFactors) {
		
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
	
	public List<int[]> subsets(List<CategoricalAttribute> byFactors) {
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
		
		int nobs = byFactors.get(0).size();
		OpenBitSet missing = new OpenBitSet(nobs);
		for (PhenotypeAttribute pa : byFactors) missing.or(pa.missing());
		
		int nfactors = byFactors.size();
		Multimap<subset, Integer> subsetMap = HashMultimap.create();
		for (int obs = 0; obs < nobs; obs++) if (!missing.fastGet(obs)) {
			int[] levels = new int[nfactors];
			int count = 0;
			for (CategoricalAttribute ca : byFactors) {
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

	public void setTraitsToTransform(List<NumericAttribute> traitsToTransform) {
		this.traitsToTransform = traitsToTransform;
	}

	public void setByFactor(List<CategoricalAttribute> byFactor) {
		this.byFactor = byFactor;
	}

	public void setLogTransform(boolean logTransform) {
		this.logTransform = logTransform;
	}

	public void setPowerTransform(boolean powerTransform) {
		this.powerTransform = powerTransform;
	}

	public void setStandardize(boolean standardize) {
		this.standardize = standardize;
	}

	public void setMyBase(BASE myBase) {
		this.myBase = myBase;
	}

	public void setPower(double power) {
		this.power = power;
	}

}
