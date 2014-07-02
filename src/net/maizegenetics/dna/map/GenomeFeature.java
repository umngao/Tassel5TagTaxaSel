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

    public static enum StrandSide {PLUS, MINUS, UNKNOWN};

    //Variables to store the information on the feature
    private String id;
    private String type;
    private int chromosome;  //Replace with a Chromosome class, or not necessary?
    private int start, stop;    //Location on the chromosome (start and stop should be inclusive)
    private StrandSide strand = StrandSide.UNKNOWN; //Strand.

    //Variables to link to parents and mychildren
    GenomeFeature parent=null;
    Multimap<String, GenomeFeature> children=null;   //Hashmap of mychildren, sorted by type (

    /**
     * Constructor to create a new GenomeFeature. Should ONLY be called by the GenomeFeatureBuilder class
     * @param myId  Unique identifier for this feature
     * @param mytype    Type of feature (gene, exon, etc)
     * @param mychr Chromosome
     * @param mystart   Start position
     * @param mystop    Stop position
     * @param mystrand  Strand (plus, minus, or unknown)
     * @param myparent  Parent feature
     * @param mychildren    Children features, in a Multimap by type
     */
    GenomeFeature(String myId, String mytype, int mychr, int mystart, int mystop, StrandSide mystrand,
                  GenomeFeature myparent, Multimap<String, GenomeFeature> mychildren){
        this.id = myId;
        this.type=mytype;
        this.chromosome=mychr;
        this.start=mystart;
        this.stop=mystop;
        this.strand=mystrand;
        this.parent=myparent;
        this.children=mychildren;
    }

    public String id(){
        return this.id;
    }

    public String type(){
        return this.type;
    }

    public int chromosome(){
        return this.chromosome;
    }

    public int start(){
        return this.start;
    }

    public int stop(){
        return this.stop;
    }

    public StrandSide strand(){
        return this.strand;
    }

    public int strandAsInt(){
        switch(strand){
            case PLUS: return 1;    //No need for break statements b/c of return
            case MINUS: return -1;
            default: return 0;
        }
    }

    public GenomeFeature parent(){
        return this.parent;
    }

    /**
     * Get all the mychildren of this feature as a single Collection
     * @return A Collection containing all mychildren GenomeFeatures of this feature
     */
    //TODO: Wrap collection into a HashSet?
    public Collection<GenomeFeature> children(){
        return children.values();
    }

    /**
     * Get a collection of all the mychildren of this feature that are annotated as a certain type
     * @param type The type of mychildren to return (gene, exon, etc)
     * @return A Collection containing all matching mychildren
     */
    //TODO: Wrap collection into a HashSet?
    public Collection<GenomeFeature> childrenOfType(String type){
        return children.get(type);
    }

    /**
     * Get the actual Multimap that stores the type->mychildren mapping for this GenomeFeature
     * @return An unmodifiable Multimap with type->child mapping (e.g., exon->GRMZ01482E02, exon->GRMZ01482E03, etc)
     */
    public Multimap<String, GenomeFeature> childrenMap(){
        return Multimaps.unmodifiableMultimap(this.children);
    }
}
