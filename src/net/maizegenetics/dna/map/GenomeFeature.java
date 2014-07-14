package net.maizegenetics.dna.map;

import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Jason Wallace on 7/2/14.
 * This class stores a generic "feature" on a genome, such as a gene, transcript, exon, miRNA, etc. The intent is for
 * this information to be read in from an external file (such as an Ensembl GFF/GTF file) and collated into a
 * GenomeFeatureMap, but it could be used in other ways.
 */
//TODO: Change int chromosome to Chromosome class? Or String to handle scaffolds?
public class GenomeFeature {

    private int start, stop;    //Location on the chromosome (start and stop should be inclusive)
    private HashMap<String, String> annotations;    //Hashmap of all annotations, stored as strings.

    GenomeFeature(HashMap<String, String> myannotations){
       this.annotations=myannotations;

       //Assign position values based on annotations. Lookup is 100-1000x faster this way than having to convert from String each time
       if(annotations.containsKey("start")){
           this.start = Integer.parseInt(annotations.get("start"));
       }
       if(annotations.containsKey("stop")){
            this.stop = Integer.parseInt(annotations.get("stop"));
       }
    }

    public String id(){
        return annotations.get("id");
    }

    public String type(){
        return annotations.get("type");
    }

    public String parentId(){
        return annotations.get("parent_id");
    }

    public String chromosome(){
        return annotations.get("chromosome");
    }

    public int start(){
        return this.start;
    }

    public int stop(){
        return this.stop;
    }

    /**
     * Returns a (shallow) copy of the Hashmap that keeps all annotations for this feature. Since the hashmap just
     * stores Strings, a shallow copy is still safe to modify b/c it won't be reflected in the original.
     * @return A copy of the Hashmap storing this feature's annotations.
     */
    public HashMap<String, String> annotations(){
        HashMap<String, String> annotationsCopy = new HashMap<>(annotations);
        return annotationsCopy;
    }

}
