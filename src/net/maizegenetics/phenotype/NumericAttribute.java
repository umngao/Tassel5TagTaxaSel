package net.maizegenetics.phenotype;

import net.maizegenetics.util.BitSet;

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
	public Object getValue(int obs) {
		return new Float(values[obs]);
	}

	@Override
	public Object getValues() {
		return values;
	}

	@Override
	public Object getSubsetOfValues(int[] obs) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public PhenotypeAttribute getSubset(int[] obs) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isMissing(int obs) {
		return missing.fastGet(obs);
	}

	@Override
	public BitSet getMissing() {
		return missing;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public int getSize() {
		return values.length;
	}

}
