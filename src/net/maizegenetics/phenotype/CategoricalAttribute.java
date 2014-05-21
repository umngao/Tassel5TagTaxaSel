package net.maizegenetics.phenotype;

import java.util.TreeSet;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import com.google.common.collect.ImmutableBiMap;

public class CategoricalAttribute implements PhenotypeAttribute {

	private final String name;
	private final int[] values;
	private final ImmutableBiMap<String, Integer> labelBimap;
	private final BitSet missing;

	public CategoricalAttribute(String name, String[] stringValues , BitSet missing) {
		this.name = name;
		int n = stringValues.length;
		this.missing = missing;
		values = new int[n];
		TreeSet<String> labelSet = new TreeSet<String>();
		for (String label : stringValues) labelSet.add(label);
		ImmutableBiMap.Builder<String, Integer> bimapBuilder = ImmutableBiMap.builder();
		int count = 0;
		for (String label : labelSet) bimapBuilder.put(label, count++);
		labelBimap = bimapBuilder.build();
	}
	
	public int getIntegerValue(int obs) {
		return values[obs];
	}
	
	public int[] getIntegerValues() {
		return values;
	}
	
	String getAttrLabel(int index) {
		return labelBimap.inverse().get(index);
	}
	
	int getIndexForAttrLabel(String label) {
		return labelBimap.get(label);
	}

	String[] getValuesAsLabels() {
		int n = values.length;
		String[] labels = new String[n];
		ImmutableBiMap<Integer, String> reverseMap = labelBimap.inverse();
		for (int i = 0; i < n; i++) {
			labels[i] = reverseMap.get(i);
		}
		return labels;
	}
	
	@Override
	public Object getValue(int obs) {
		return labelBimap.inverse().get(values[obs]);
	}

	@Override
	public Object getValues() {
		return getValuesAsLabels();
	}

	@Override
	public Object getSubsetOfValues(int[] obs) {
		int n = obs.length;
		String[] labels = new String[n];
		ImmutableBiMap<Integer, String> reverseMap = labelBimap.inverse();
		for (int i = 0; i < n; i++) {
			labels[i] = reverseMap.get(obs[i]);
		}
		return labels;
	}

	@Override
	public PhenotypeAttribute getSubset(int[] obs) {
		int n = obs.length;
		String[] labels = new String[n];
		OpenBitSet subMissing = new OpenBitSet(n);
		ImmutableBiMap<Integer, String> reverseMap = labelBimap.inverse();
		for (int i = 0; i < n; i++) {
			labels[i] = reverseMap.get(obs[i]);
			if (missing.fastGet(i)) subMissing.fastSet(i);
		}
		return new CategoricalAttribute(name, labels, subMissing);
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
