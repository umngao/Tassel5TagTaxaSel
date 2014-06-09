package net.maizegenetics.util;

import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.SetMultimap;

import java.util.Map;

/**
 * Not sure we want to do this, but if we do then we need to use the myAnno approach, as it massively reduces the overhead
 * for small number of annotations.
 *
 * @author Ed Buckler
 */
public class AbstractAnnotation implements GeneralAnnotation {
    public static final Object[] EMPTYOBJ=new Object[0];
    public static final String[] EMPTYSTR=new String[0];
    public static final double[] EMPTYdouble=new double[0];
    //Custom annotation are stored in the map
    private final ImmutableMultimap<String, Object> myAnnoMap;

  //  private final Map.Entry<String, String>[] myAnno;


    public AbstractAnnotation(ImmutableMultimap<String, Object> myAnnoMap) {
        this.myAnnoMap=myAnnoMap;
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        try{return myAnnoMap.get(annoName).toArray();}
        catch(Exception e) {return EMPTYOBJ;}
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        try{return myAnnoMap.get(annoName).toArray(new String[0]);}
        catch(Exception e) {return EMPTYSTR;}
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        try{
            Object[] o=myAnnoMap.get(annoName).toArray();
            if((o == null)||(!(o[0] instanceof Number))) return EMPTYdouble;
            double[] d=new double[o.length];
            int i=0;
            for (Object o1 : o) {d[i++]=((Number)o1).doubleValue();}
            return d;
        }catch(Exception e) {
            return EMPTYdouble;
        }
    }

    @Override
    public String getConsensusAnnotation(String annoName) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Map.Entry[] getAllAnnotationEntries() {
        return new Map.Entry[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public SetMultimap<String, String> getAnnotationAsMap() {
        ImmutableSetMultimap.Builder<String,String> result=new ImmutableSetMultimap.Builder<String,String>()
                .orderKeysBy(Ordering.natural()).orderValuesBy(Ordering.natural());
        for (Map.Entry<String, Object> en : getAllAnnotationEntries()) {
            result.put(en.getKey(),en.getValue().toString());
        }
        return result.build();
    }

    @Override
    public boolean isAnnotatedWithValue(String annoName, String annoValue) {
        return false;
    }
}
