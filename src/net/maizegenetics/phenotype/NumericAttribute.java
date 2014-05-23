package net.maizegenetics.phenotype;



import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class NumericAttribute implements PhenotypeAttribute {
	private final String name;
	private final float[] values;
	private final BitSet missing;

	public NumericAttribute(String name, float[] values, BitSet missing) {
		this.name = name;
		this.values = values;
		this.missing = missing;
	}
	
	public float getFloatValue(int obs) {
		return values[obs];
	}
	
	public float[] getFloatValues() {
		return values;
	}
	
	@Override
	public Object value(int obs) {
		return new Float(values[obs]);
	}

	@Override
	public Object allValues() {
		return values;
	}

	@Override
	public PhenotypeAttribute subset(int[] obs) {
		int n = obs.length;
		float[] valueSubset = new float[n];
		OpenBitSet missingSubset = new OpenBitSet(n);
		for (int i = 0; i < n; i++) {
			valueSubset[i] = values[obs[i]];
			if (missing.fastGet(obs[i])) missingSubset.fastSet(i);
		}
		return new NumericAttribute(name, valueSubset, missingSubset);
	}

	@Override
	public boolean isMissing(int obs) {
		return missing.fastGet(obs);
	}

	@Override
	public BitSet missing() {
		return missing;
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public int size() {
		return values.length;
	}

}
