/*
 *  GeneralAnnotationStorage
 * 
 *  Created on Aug 5, 2014
 */
package net.maizegenetics.util;

import ch.systemsx.cisd.hdf5.IHDF5Reader;

import com.google.common.collect.SetMultimap;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class GeneralAnnotationStorage implements GeneralAnnotation {

    private static final double[] EMPTY_DOUBLE_ARRAY = new double[0];

    private static final int MAX_CACHE_SIZE = 1_000_000;

    private static final Map<Map.Entry<String, String>, Map.Entry<String, String>> CACHE = Collections.synchronizedMap(new LinkedHashMap<Map.Entry<String, String>, Map.Entry<String, String>>((3 * MAX_CACHE_SIZE) / 2) {

        @Override
        protected boolean removeEldestEntry(Map.Entry<Map.Entry<String, String>, Map.Entry<String, String>> eldest) {
            return size() > MAX_CACHE_SIZE;
        }

    });

    static Map.Entry<String, String> getCanonicalAnnotation(String key, String value) {
        Map.Entry<String, String> temp = new AbstractMap.SimpleImmutableEntry<>(key, value);
        Map.Entry<String, String> entry = CACHE.putIfAbsent(temp, temp);
        return (entry == null) ? temp : entry;
    }

    private final Map.Entry<String, String>[] myAnnotations;

    private GeneralAnnotationStorage(Builder builder) {
        myAnnotations = (Map.Entry<String, String>[]) new Map.Entry<?, ?>[builder.myAnnotations.size()];
        for (int i = 0; i < builder.myAnnotations.size(); i++) {
            myAnnotations[i] = builder.myAnnotations.get(i);
        }
    }

    public static Builder getBuilder() {
        return new Builder();
    }

    public static GeneralAnnotationStorage getFromHDF5(IHDF5Reader reader) {
        return null;
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        List<Object> result = new ArrayList<>(1);
        for (Map.Entry<String, String> me : myAnnotations) {
            if (me.getKey().equals(annoName)) {
                result.add(me.getValue());
            }
        }
        return result.toArray();
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        List<String> result = new ArrayList<>(1);
        for (Map.Entry<String, String> me : myAnnotations) {
            if (me.getKey().equals(annoName)) {
                result.add(me.getValue());
            }
        }
        return result.toArray(new String[result.size()]);
    }

    @Override
    public String getConsensusAnnotation(String annoName) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        try {
            ArrayList<Double> result = new ArrayList<>(1);
            for (Map.Entry<String, String> me : myAnnotations) {
                if (me.getKey().equals(annoName)) {
                    result.add(Double.parseDouble(me.getValue()));
                }
            }
            if (result.isEmpty()) {
                return EMPTY_DOUBLE_ARRAY;
            }
            double[] d = new double[result.size()];
            for (int i = 0; i < result.size(); i++) {
                d[i] = result.get(i);
            }
            return d;
        } catch (Exception e) {
            return EMPTY_DOUBLE_ARRAY;
        }
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Map.Entry<String, String>[] getAllAnnotationEntries() {
        return Arrays.copyOf(myAnnotations, myAnnotations.length);
    }

    @Override
    public SetMultimap<String, String> getAnnotationAsMap() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isAnnotatedWithValue(String annoName, String annoValue) {
        for (Map.Entry<String, String> me : myAnnotations) {
            if (me.getKey().equals(annoName) && me.getValue().equals(annoValue)) {
                return true;
            }
        }
        return false;
    }

    public static class Builder {

        private final List<Map.Entry<String, String>> myAnnotations = new ArrayList<>(0);

        private Builder() {
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnnotation(String key, String value) {
            Map.Entry<String, String> ent = getCanonicalAnnotation(key, value);
            myAnnotations.add(ent);
            return this;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnnotation(String key, Number value) {
            Map.Entry<String, String> ent = getCanonicalAnnotation(key, value.toString());
            myAnnotations.add(ent);
            return this;
        }

        public GeneralAnnotationStorage build() {
            if (myAnnotations.size() == 0) {
                return null;
            }
            Collections.sort(myAnnotations, new Comparator<Map.Entry<String, String>>() {
                public int compare(Map.Entry<String, String> s1, Map.Entry<String, String> s2) {
                    int keyComp = s1.getKey().compareTo(s2.getKey());
                    if (keyComp != 0) {
                        return keyComp;
                    }
                    return s1.getValue().compareTo(s2.getValue());
                }
            });
            return new GeneralAnnotationStorage(this);
        }
    }
}
