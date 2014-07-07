package net.maizegenetics.dna.map;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.map.GenomeFeature.StrandSide;

/**
 * Created by jgw87 on 7/2/14.
 * Builder class to create a GenomeFeature.
 * Note: There are no methods add children without an associated type (exon, gene, etc). Please do not create them.
 */
//TODO: There's a flaw in this, in that to link to parents and children they must already be built, but when constructing
    //a map that's not possible. Switch to String names and a link to the master lookup table instead?
public class GenomeFeatureBuilder {

    //Variables to store the information on the feature
    private String myId=null;
    private String mytype =null;
    private int mychromosome =-1;  //Replace with a Chromosome class, or not necessary?
    private int mystart =-1, mystop =-1;    //Location on the mychromosome (mystart and mystop should be inclusive)
    private StrandSide mystrand =StrandSide.UNKNOWN; //Strand.

    //Variables to link to parents and mychildren
    GenomeFeature myparent =null;
    Multimap<String, GenomeFeature> mychildren = HashMultimap.create();

    /**
     * Generic constructor which does nothing special
     */
    public GenomeFeatureBuilder(){}

    /**
     * Constructor to build a new feature off of an existing one
     * @param feature
     */
    public GenomeFeatureBuilder(GenomeFeature feature){
        this.mytype =feature.type();
        this.mychromosome =feature.chromosome();
        this.mystart =feature.start();
        this.mystop =feature.stop();
        this.mystrand =feature.strand();
        this.myparent =feature.parent();

        //Since the multimap returned by GenomeFeature.childrenmap() is immutable, copy it into a new multimap
        this.mychildren = HashMultimap.create();
        Multimap <String, GenomeFeature> childrenMap = feature.childrenMap();
        for(String type: childrenMap.keySet()){
            for(GenomeFeature child: childrenMap.get(type)){
                mychildren.put(type, child);
            }
        }
    }

    /**
     * Public accessor method to get a new GenomeFeatureBuilder based off an existing GenomeFeature
     * TODO: Get this matching other getInstance methods better; doesn't seem to fit, so commented out for now
     * //@param feature
     * //@return
     */
    /*public static GenomeFeatureBuilder getInstance(GenomeFeature feature){
        return new GenomeFeatureBuilder(feature);
    }*/

    public GenomeFeature build(){
        return new GenomeFeature(myId, mytype, mychromosome, mystart, mystop, mystrand, myparent, mychildren);
    }

    public GenomeFeatureBuilder id(String id){
        myId=id;
        return this;
    }

    public GenomeFeatureBuilder type(String type){
        mytype=type;
        return this;
    }

    public GenomeFeatureBuilder chromosome(int chromosome){
        mychromosome=chromosome;
        return this;
    }

    public GenomeFeatureBuilder start(int start){
        mystart=start;
        return this;
    }

    public GenomeFeatureBuilder stop(int stop){
        mystop=stop;
        return this;
    }

    public GenomeFeatureBuilder strand(StrandSide strand){
        mystrand=strand;
        return this;
    }

    /**
     * Assign strand based on an interger value. 1 is plus strand, -1 is minus strand. Anything else is unknown.
     * @param strand Strand side (1 or -1)
     * @return This builder
     */
    public GenomeFeatureBuilder strand(int strand){
        switch(strand){
            case 1: mystrand=StrandSide.PLUS; break;
            case -1: mystrand=StrandSide.MINUS; break;
            default: mystrand=StrandSide.UNKNOWN; break;
        }
        return this;
    }

    public GenomeFeatureBuilder parent(GenomeFeature parent){
        myparent=parent;
        return this;
    }

    /**
     * Add a child. Since this
     * @param type  The type of this child (gene, exon, etc)
     * @param child The child GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureBuilder child(String type, GenomeFeature child){
        mychildren.put(type, child);
        return this;
    }

    /**
     * Add a group of children, all of the same type. The children must be in some sort of structure that implements
     * the Iterable interface (usually some sort of Collection, such as a List, Set, etc.)
     * @param type The type of all children to be added (gene, exon, etc)
     * @param children An Iterable-compatible collection of the children to be added
     * @return This builder
     */
    public GenomeFeatureBuilder children(String type, Iterable<GenomeFeature> children){
        for(GenomeFeature c: children){
            mychildren.put(type, c);
        }
        return this;
    }
}
