package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import com.google.common.collect.ImmutableBiMap;

public class CategoricalAttribute implements PhenotypeAttribute {

	private final String name;
	private final int[] values;
	private final ImmutableBiMap<String, Integer> labelBimap;
	private final BitSet missing;
	
	private static final List<ATTRIBUTE_TYPE> myAllowedTypes;
	static{
		myAllowedTypes = new ArrayList<ATTRIBUTE_TYPE>();
		myAllowedTypes.add(ATTRIBUTE_TYPE.factor);
	}

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
	
	public int intValue(int obs) {
		return values[obs];
	}
	
	public int[] allIntValues() {
		return values;
	}
	
	public String getAttrLabel(int index) {
		return labelBimap.inverse().get(index);
	}
	
	public int indexForAttrLabel(String label) {
		return labelBimap.get(label);
	}

	public String label(int obs) {
		return labelBimap.inverse().get(obs);
	}
	
	public String[] allLabels() {
		int n = values.length;
		String[] labels = new String[n];
		ImmutableBiMap<Integer, String> reverseMap = labelBimap.inverse();
		for (int i = 0; i < n; i++) {
			labels[i] = reverseMap.get(i);
		}
		return labels;
	}
	
	public List<String> labelList() {
		int n = labelBimap.size();
		ArrayList<String> labelList = new ArrayList<>();
		ImmutableBiMap<Integer, String> reverseMap = labelBimap.inverse();
		for (int i = 0; i < n; i++) labelList.add(reverseMap.get(i));
		return labelList;
	}
	
	@Override
	public Object value(int obs) {
		return labelBimap.inverse().get(values[obs]);
	}

	@Override
	public Object allValues() {
		return allLabels();
	}

	@Override
	public PhenotypeAttribute subset(int[] obs) {
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

	@Override
	public List<ATTRIBUTE_TYPE> getCompatibleTypes() {
		return myAllowedTypes;
	}

	@Override
	public boolean isTypeCompatible(ATTRIBUTE_TYPE type) {
		return myAllowedTypes.contains(type);
	}

}
