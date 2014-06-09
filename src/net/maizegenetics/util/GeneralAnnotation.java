package net.maizegenetics.util;

import com.google.common.collect.SetMultimap;

import java.util.Map;

/**
 * Provide generalized annotations (descriptors) for a taxon or site.
 *
 */
public interface GeneralAnnotation {

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public Object[] getAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public String[] getTextAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public String getConsensusAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new double[0])
     */
    public double[] getQuantAnnotation(String annoName);

    /**
     * Returns average annotation for a given annotation key
     * @param annoName annotation key
     * @return average value (if not present - return Double.NaN)
     */
    public double getAverageAnnotation(String annoName);

    /**
     * Returns all annotation Map.Entries.
     * @return array of Map.Entry
     */
    public Map.Entry<String, String>[] getAllAnnotationEntries();

    /**
     * Returns all annotations in a TreeMap.
     * @return Map of annotations
     */
    public SetMultimap<String, String> getAnnotationAsMap();

    //should we provide methods, to average the quantitative annotations, the first annotation
    //

    /**
     * Returns whether the entity contains the annotation with the specified value.  If either the annotation or the
     * value is missing false is return
     * @param annoName annotation key
     * @param annoValue annotation value;
     * @return
     */
    public boolean isAnnotatedWithValue(String annoName, String annoValue);

}
