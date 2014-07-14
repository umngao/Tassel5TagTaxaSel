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

   // public static enum StrandSide {PLUS, MINUS, UNKNOWN}; //Strand now a miscellaneous feature

    //Variables to store the position information on the feature
    //private String id;
    //private String type;
    //private String parentId;
    //private Chromosome chromosome;  //Replace with a Chromosome class, or not necessary?
    //private String chromosome;
    private int start, stop;    //Location on the chromosome (start and stop should be inclusive)
    private HashMap<String, String> annotations;    //Hashmap of all annotations, stored as strings.
    //private StrandSide strand = StrandSide.UNKNOWN; //Strand.

    //TODO: Do some speed testing to see how much difference it makes to store start-stop in own variables instead of just converting the Hashmap values each time

    /*//Variables to link to parents and mychildren
    DEPRECATED - GenomeFeatureMap uses an explicit graph structure instead
    GenomeFeature parent=null;
    Multimap<String, GenomeFeature> children=null;   //Hashmap of mychildren, sorted by type (*/

    /**
     * Constructor to create a new GenomeFeature. Should ONLY be called by the GenomeFeatureBuilder class
    // * @param myId  Unique identifier for this feature
    // * @param mytype    Type of feature (gene, exon, etc)
    // * @param mychr Chromosome
    // * @param mystart   Start position
    // * @param mystop    Stop position
     * @param myannotations  Hashmap of _all_ annotations, stored as Strings
     * //@param mystrand  Strand (plus, minus, or unknown)
     * //@param myparent  Parent feature
     * //@param mychildren    Children features, in a Multimap by type
     */
    /*GenomeFeature(String myId, String mytype, Chromosome mychr, int mystart, int mystop, HashMap<String, String> myannotations, String myParentId){
        this.id = myId;
        this.type=mytype;
        this.chromosome=mychr;
        this.start=mystart;
        this.stop=mystop;
        //this.strand=mystrand;
        this.annotations = myannotations;
        this.parentId=myParentId;
        //this.parent=myparent;
        //this.children=mychildren;
    }*/

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

    /*public Chromosome chromosome(){
        return this.chromosome;
    }*/

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

    /*public StrandSide strand(){   //Deprecated; strand isn't important for our things. Keep as a misc. annotation
        return this.strand;
    }*/

    /*public int strandAsInt(){
        switch(strand){
            case PLUS: return 1;    //No need for break statements b/c of return
            case MINUS: return -1;
            default: return 0;
        }
    }*/

    //ALL PARENT-CHILDREN FEATURES DEPRECATED

    /*public GenomeFeature parent(){
        return this.parent;
    }*/

    /**
     * Get all the mychildren of this feature as a single Collection
     * @return A Collection containing all mychildren GenomeFeatures of this feature
     */
    //TODO: Wrap collection into a HashSet?
    /*public Collection<GenomeFeature> children(){
        return children.values();
    }*/

    /**
     * Get a collection of all the mychildren of this feature that are annotated as a certain type
     * @param type The type of mychildren to return (gene, exon, etc)
     * @return A Collection containing all matching mychildren
     */
    //TODO: Wrap collection into a HashSet?
    /*public Collection<GenomeFeature> childrenOfType(String type){
        return children.get(type);
    }*/

    /**
     * Get the actual Multimap that stores the type->mychildren mapping for this GenomeFeature
     * @return An unmodifiable Multimap with type->child mapping (e.g., exon->GRMZ01482E02, exon->GRMZ01482E03, etc)
     */
    /*public Multimap<String, GenomeFeature> childrenMap(){
        return Multimaps.unmodifiableMultimap(this.children);
    }*/
}
