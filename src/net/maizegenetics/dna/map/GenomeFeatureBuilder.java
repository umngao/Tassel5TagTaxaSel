package net.maizegenetics.dna.map;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.SetMultimap;
//import net.maizegenetics.dna.map.GenomeFeature.StrandSide;
import org.apache.commons.math.exception.NumberIsTooSmallException;
import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by jgw87 on 7/2/14.
 * Builder class to create a GenomeFeature.
 * Note: There are no methods add children without an associated type (exon, gene, etc). Please do not create them.
 */
//TODO: Currently somewhat messy back-and-forth conversion for start and stop (can go String-int-String-int). Find a better way?
//TODO: Set start & stop to 0 initially, not -1? (Implicitly required now)
//TODO: Make id, chrom, start, stop, & type all required data? Makes sense given the uses of them.
public class GenomeFeatureBuilder {

    private static final Logger myLogger = Logger.getLogger(GenomeFeatureBuilder.class);

    //Variables to store the information on the feature
    //private String myId=null;
    //private String mytype =null;
    //private String myParentId=null;
    //private Chromosome mychromosome = null;
    //private String mychromosome = null;
    private int mystart =-1, mystop =-1;    //Location on the mychromosome (mystart and mystop should be inclusive). Negative by default to flag as unassigned
    private HashMap<String, String> myannotations = null;
    //private StrandSide mystrand =StrandSide.UNKNOWN; //Strand.

    //Variables to link to parents and mychildren - DEPRECATED. GenomeFeatureMapBuilder uses explicit graph instead
    //GenomeFeature myparent =null;
    //Multimap<String, GenomeFeature> mychildren = HashMultimap.create();

    /**
     * Generic constructor which does nothing special
     */
    public GenomeFeatureBuilder(){
        myannotations = new HashMap<String, String>();
    }

    /**
     * Constructor to build a new feature off of an existing one
     * @param feature
     */
    public GenomeFeatureBuilder(GenomeFeature feature){
        //this.mytype =feature.type();
        //this.mychromosome =feature.chromosome();
        this.mystart =feature.start();
        this.mystop =feature.stop();
        //this.mystrand =feature.strand();
        this.myannotations = feature.annotations();
        //this.myParentId= feature.parentId();
        //this.myparent =feature.parent();

        //Since the multimap returned by GenomeFeature.childrenmap() is immutable, copy it into a new multimap
        /*this.mychildren = HashMultimap.create();
        Multimap <String, GenomeFeature> childrenMap = feature.childrenMap();
        for(String type: childrenMap.keySet()){
            for(GenomeFeature child: childrenMap.get(type)){
                mychildren.put(type, child);
            }
        }*/
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
        validateData();
        //return new GenomeFeature(myId, mytype, mychromosome, mystart, mystop, myannotations, myParentId);
        return new GenomeFeature(myannotations);
    }

    private void validateData(){
        //System.out.println("Building...\n\tStart=" + mystart + "_" + myannotations.get("start") + "\n\tStop=" + mystop + "_" + myannotations.get("stop"));

        //Test that feature has a personal ID
        if(!myannotations.containsKey("id")){
            throw new UnsupportedOperationException("GenomeFeatureBuilder: Cannot build a feature without a personal identifier (field 'id')");
        }

        //Test if start or stop is negative (eg, if was never assigned)
        if(mystart < 0){
            throw new UnsupportedOperationException("GenomeFeatureBuilder: Start coordinate is negative for " +
                    myannotations.get("id") + ": " + mystart + " (possibly unassigned?)");
        }
        if(mystop < 0){
            throw new UnsupportedOperationException("GenomeFeatureBuilder: Stop coordinate is negative for " +
                    myannotations.get("id") + ": " + mystop + " (possibly unassigned?)");
        }

        //Test that start is less than stop
        if(mystart > mystop){
            throw new UnsupportedOperationException("GenomeFeatureBuilder: Start coordinate is greater than stop " +
                    "coordinate for " + myannotations.get("id") + ": " +  mystart + " vs " + mystop);
        }
    }

    public GenomeFeatureBuilder id(String id){
        return addAnnotation("id", id);
    }

    public GenomeFeatureBuilder type(String type){
        return addAnnotation("type", type);
    }

    public GenomeFeatureBuilder parentId(String parentId){
        return addAnnotation("parent_id", parentId);
    }

    public GenomeFeatureBuilder chromosome(Chromosome chr){
        return addAnnotation("chromosome", chr.getName());
        //mychromosome=chr;
        //return this;
    }

    public GenomeFeatureBuilder chromosome(String chr){
        return addAnnotation("chromosome", chr);
        //mychromosome=new Chromosome("" + chr);
        //return this;
    }

    public GenomeFeatureBuilder chromosome(int chr){
        return addAnnotation("chromosome", "" + chr);
    }

    public GenomeFeatureBuilder start(int start){
        if(start >=0){
            return addAnnotation("start", "" + start);
        }else{
            throw new NumberFormatException("Start position must be greater than zero. Got " + start);
        }
    }

    public GenomeFeatureBuilder start(String start){
        return this.start(Integer.parseInt(start));
    }

    public GenomeFeatureBuilder stop(int stop){
        if(stop >=0){
            return addAnnotation("stop", "" + stop);
        }else{
            throw new NumberFormatException("Stop position must be greater than zero. Got " + stop);
        }
        //return this;
    }

    public GenomeFeatureBuilder stop(String stop){
        return this.stop(Integer.parseInt(stop));
    }


    public GenomeFeatureBuilder addAnnotation(String key, String value){
        key = synonymizeKeys(key);
        myannotations.put(key, value);  //All annotations kept in the hash

        //System.out.println("Adding "+ key + " with value " + value);

        //Specific handling for start-stop positions because are used in validation at build time
        switch(key){
            case "start": mystart = Integer.parseInt(value); break;
            case "stop": mystop = Integer.parseInt(value); break;
        }

        return this;
    }

    /**
     * Method that takes common synonyms of annotation types and standardizes them according to the following rules:
     * (1) Make lowercase
     * (2) Standardize according to following rules. (Any not on this list are returned as just lowercased)
     *    name, id -> id
     *    chr, chrom, chromosome -> chromosome
     *    stop, end -> stop
     *    parentid, parent_id, parent -> parent_id
     * @param key
     * @return
     */
    public String synonymizeKeys(String key){
        key.toLowerCase(Locale.ENGLISH);
        switch(key){
            case "name":
            case "id" : return "id";

            case "chr":
            case "chrom":
            case "chromosome": return "chromosome";

            case "end":
            case "stop": return "stop";

            case "parent":
            case "parentid":
            case "parent_id": return "parent_id";

            default: return key;
        }
    }

    /*public GenomeFeatureBuilder strand(StrandSide strand){
        mystrand=strand;
        return this;
    }*/

    /**
     * Assign strand based on an integer value. 1 is plus strand, -1 is minus strand. Anything else is unknown.
     * @param strand Strand side (1 or -1)
     * @return This builder
     */
    /*public GenomeFeatureBuilder strand(int strand){
        switch(strand){
            case 1: mystrand=StrandSide.PLUS; break;
            case -1: mystrand=StrandSide.MINUS; break;
            default: mystrand=StrandSide.UNKNOWN; break;
        }
        return this;
    }*/

    /**
     * Assign strand based on an character value. '+' for plus, '-' for minus, anything else is unknown
     * @param strand Strand side (1 or -1)
     * @return This builder
     */
    /*public GenomeFeatureBuilder strand(char strand){
        switch(strand){
            case '+': mystrand=StrandSide.PLUS; break;
            case '-': mystrand=StrandSide.MINUS; break;
            default: mystrand=StrandSide.UNKNOWN; break;
        }
        return this;
    }*/

    /**
     * Assign strand based on the first character in a string. '+' for plus, '-' for minus, anything else is unknown.
     * Note that only the first character is tested; all others are ignored.
     * @param strand Strand side (1 or -1)
     * @return This builder
     */
    /*public GenomeFeatureBuilder strand(String strand){
        return this.strand(strand.charAt(0));
    }*/


    /**
     * Create a GenomeFeature from a line of GFF file. This method is modified from the BioJava source code for the same
     * purpose, in biojava3-genome/src/main/java/org/biojava3/genome/GFF3Reader.java
     * @param line A single line from a GFF file as a string
     */
    public GenomeFeatureBuilder parseGffLine(String line){
        //Field identifiers for GFF format
        int gffSeq=0, gffSource=1, gffFeatureType=2, gffStart=3, gffStop=4, gffScore=5, gffStrand=6, gffFrame=7, gffAttributes=8;

        //Get all the easy data stored in its own fields
        String[] tokens = line.split("\t");
        this.chromosome(tokens[gffSeq].trim());
        this.type(tokens[gffFeatureType].trim());
        this.start(tokens[gffStart]);
        this.stop(tokens[gffStop]);
        addAnnotation("strand", tokens[gffStrand].trim());

        //Extract the parent from the attributes field
        String parentID = getParentFromGffAttributes(tokens[gffAttributes]);
        this.parentId(parentID);

        //Extract the unique identifier for this feature. If none, build one from available info
        String myID = getFeatureIdFromGffAttributes(tokens[gffAttributes]);
        if(myID == null){
            myID = myannotations.get("type") + "_" + myannotations.get("chromosome") + "_" + mystart + "_" + mystop;
        }
        this.id(myID);

        return this;
    }

    /**
     * Parse a GFF attribute field to identify the parent of the current GenomeFeature. Tricky b/c of the different ways it
     * can be represented. There's a hierarchy of accepted answers, with 'parent_id', 'Parent=', 'transcript_id', and
     * 'gene_id' taken in that order. If nothing is found, returns an empty string ("")
     * @param attributes The string from the attribute field of the GFF file
     * @return The parent's ID string
     */
    public static String getParentFromGffAttributes(String attributes){
        //Match the appropriate string with regular expressions
        Matcher matcher;

        //Pattern for 'parent_id "GRMZM2G005232"'
        matcher = Pattern.compile("parent_id \"(\\w+)\"").matcher(attributes);
        if(matcher.find()){
            return matcher.group(1);
        }

        //Pattern for 'Parent=GRMZM2G005232' and  'Parent=gene:GRMZM2G005232'
        matcher = Pattern.compile("Parent=(\\w+:){0,1}(\\w+)").matcher(attributes);
        if(matcher.find()){
            return matcher.group(2);
        }

        //Pattern for 'gene_id "GRMZM2G005232"'
        matcher = Pattern.compile("transcript_id \"(\\w+)\"").matcher(attributes);
        if(matcher.find()){
            return matcher.group(1);
        }

        //Pattern for 'transcript_id "GRMZM2G005232_T01"'
        matcher = Pattern.compile("gene_id \"(\\w+)\"").matcher(attributes);
        if(matcher.find()){
            return matcher.group(1);
        }

        return "";
    }

    /**
     * Parse a GFF attribute field to identify the name of the current GenomeFeature. Looks for 'ID=' and 'Name=' fields
     * @param attributes The string from the attribute field of the GFF file. REturns null if not found
     * @return The feature's ID string
     */
    public static String getFeatureIdFromGffAttributes(String attributes) {
        //Match the appropriate string with regular expressions
        Matcher matcher;

        //Pattern for 'ID=GRMZM2G005232' and  'ID=gene:GRMZM2G005232'
        matcher = Pattern.compile("(Name|ID)=(\\w+:){0,1}(\\w+)").matcher(attributes);
        if(matcher.find()){
            return matcher.group(3);
        }

        return null;
    }

    /*public GenomeFeatureBuilder parent(GenomeFeature parent){
        myparent=parent;
        return this;
    }*/

    /**
     * Add a child. Since this
     * @param type  The type of this child (gene, exon, etc)
     * @param child The child GenomeFeature to be added
     * @return This builder
     */
    /*public GenomeFeatureBuilder child(String type, GenomeFeature child){
        mychildren.put(type, child);
        return this;
    }*/

    /**
     * Add a group of children, all of the same type. The children must be in some sort of structure that implements
     * the Iterable interface (usually some sort of Collection, such as a List, Set, etc.)
     * @param type The type of all children to be added (gene, exon, etc)
     * @param children An Iterable-compatible collection of the children to be added
     * @return This builder
     */
   /*public GenomeFeatureBuilder children(String type, Iterable<GenomeFeature> children){
        for(GenomeFeature c: children){
            mychildren.put(type, c);
        }
        return this;
    }*/
}
