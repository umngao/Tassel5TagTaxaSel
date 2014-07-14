package net.maizegenetics.dna.map;

import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import net.maizegenetics.util.DirectedGraph;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Created by jgw87 on 7/2/14.
 * A Builder class to create a GenomeFeatureMap to identify genomic features. Can build piecemeal or read in from a file
 * For now, not implementing functionality to make a new builder from an existing map.
 */
//TODO: Add functionality for reading in a precompiled map
public class GenomeFeatureMapBuilder {

    //Root of the GenomeFeature tree - DEPRECATED IN FAVOR OF AN ACTUAL GRAPH
    //GenomeFeature root = null;

    private static final Logger myLogger = Logger.getLogger(GenomeFeatureMapBuilder.class);

    //Graph of all the genome features connecting to each other
    DirectedGraph<GenomeFeature> featureTree = null;

    //Lookups to identify GenomeFeatures by their name, location, and type (exon, gene, etc)
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup = null; // Complex.
    private Multimap<String, GenomeFeature> typeLookup = null;

    //Helper variables used to store information to build the feature map
    HashSet<String> chromosomes = new HashSet<>();



    public GenomeFeatureMap build(){
        //buildGenomeTree();
        buildLocationLookup();
        return new GenomeFeatureMap(nameLookup, typeLookup, locationLookup, featureTree);
    }

    private void buildGenomeTree(){
        //TODO: Implement

        //TODO: Build chromosomes first

        /*//Pull out the "genome" feature to serve as root. If it doesn't exist, create it.
        if(typeLookup.keySet().contains("genome")){
            Collection<GenomeFeature> mygenome = typeLookup.get("genome");
            if(mygenome.size() >1){ //If more than one feature annotated as "genome", throw an error
                StringBuilder types = new StringBuilder();
                for(GenomeFeature f:mygenome){
                    types.append("\n\t" + f.id());
                }
                throw new UnsupportedOperationException("Error: Attempt to build a GenomeFeatureMap with more than one feature annotated as type 'genome':" + types);
            }
            root = mygenome.toArray(new GenomeFeature[0])[0];
        }else{
            GenomeFeatureBuilder rootBuilder = new GenomeFeatureBuilder()
                    .parent(null).chromosome(-1).id("GENOME").type("genome");

            root = new GenomeFeatureBuilder().parent(null).build();
        }*/
    }

    private void buildLocationLookup(){
        //Initialize each chromosome as a new Rangemap
        for(String c: chromosomes){
            RangeMap<Integer, HashSet<GenomeFeature>> newmap = TreeRangeMap.create();
            locationLookup.put(c, newmap);
        }

        for(GenomeFeature feature: nameLookup.values()){
            RangeMap mymap = locationLookup.get(feature.chromosome());
            addFeatureToRangemap(mymap, feature);
            //TODO: I think the logic here is faulty
        }
    }

    /**
     * Add a GenomeFeature to a RangeMap, stacking it on top of any existing features instead of overwriting them
     * (how RangeMap behaves natively)
     * @param masterMap Rangemap object to be modified
     * @param feature GenomeFeature to be added. Only start and stop position are checked
     * @return
     */
    public static void addFeatureToRangemap(
            RangeMap<Integer, HashSet<GenomeFeature>>  masterMap, GenomeFeature feature){

        Range<Integer> span = Range.closed(feature.start(), feature.stop());

        //First, extract out subrange (to preserve annotations already there)
        RangeMap subrange = masterMap.subRangeMap(span);

        //Next, set entire range equal to the new feature (so that any previously empty regions now map to it)
        HashSet<GenomeFeature> newset = new HashSet<>();
        newset.add(feature);
        masterMap.put(span, newset);

        //Take the extracted, existing annotations and add back in, adding this feature to the list
        Map<Range<Integer>, HashSet<GenomeFeature>> subranges = subrange.asMapOfRanges();
        for(Range r: subranges.keySet()){
            HashSet<GenomeFeature> tempset = subranges.get(r);
            tempset.add(feature);
            masterMap.put(r, tempset);
        }

        //return masterMap;
    }

    /**
     * Adds a GenomeFeature to the map. If the feature's unique ID has already been loaded, it throws an UnsupportedOperationException
     * @param feature The GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureMapBuilder addFeature(GenomeFeature feature){
        String ID=feature.id();
        if(nameLookup.containsKey(ID)){
            throw new UnsupportedOperationException("Error: Attempt to add a GenomeFeature whose unique ID is already loaded: "
            + ID);
        }
        return this.addOrReplaceFeature(feature);   //Done to avoid code duplication
    }

    /**
     * Replaces the GenomeFeature with a specified ID with a new one. If unique ID isn't already in the map, it throws
     * an UnsupportedOperationException.
     * @param feature The GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureMapBuilder replaceFeature(GenomeFeature feature){
        String ID=feature.id();
        if(!nameLookup.containsKey(ID)){
            throw new UnsupportedOperationException("Error: Attempt to replace a GenomeFeature whose unique ID has not been loaded yet: "
                    + ID);
        }
        return this.addOrReplaceFeature(feature);   //Done to avoid code duplication
    }

    /**
     * Adds a GenomeFeature to the map, regardless of whether or not it's already been added. This method throws no
     * warnings if you'll overwrite existing data, so use it with caution.
     * @param feature The GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureMapBuilder addOrReplaceFeature(GenomeFeature feature){
        nameLookup.put(feature.id(), feature);
        typeLookup.put(feature.type(), feature);
        chromosomes.add(feature.chromosome());
        return this;
    }

    /**
     * Load in data from a GFF (Gene Feature Format) file. Since GFF files have only a loose standard for the 9th column,
     * this involves several ad-hoc heuristics about what things to look for. As such, it is not the preferred way to
     * read in annotations. (That is JSON or tab-delimited format.)
     *
     * This method does not build the map, so you can string
     * multiple calls together (if, for example, you have different annotations in different files)
     *
     *
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromGtfFile(String filename){
        myLogger.warn("GenomeFeatureMapBuilder - Loading genome annotations from GFF file. Will try to parse annotations " +
                "field as best as possible. (JSON or tab-delimited formats are preferred.)");
        try {
            BufferedReader reader = Utils.getBufferedReader(filename);
            String line=reader.readLine();
            while(line != null ){
                if(line.startsWith("#")){    //Skip comment lines
                    line=reader.readLine();
                    continue;
                }
                GenomeFeature newFeature = new GenomeFeatureBuilder().parseGffLine(line).build();
                addFeature(newFeature);
                line=reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return this;
    }


    /**
     * Load in data from a JSON-formatted file. JSON format is defined at http://www.json.org/, and consists of structured
     * key-value pairs. For genome features, the key is the name of an attribute and the value is (obviously) its value.
     * (For example: "chromosome":1). Common attributes (keys) are listed below. Although only the "id" attribute is
     * required, a feature is pretty useless without some sort of positional information (chromosome, start/stop, etc.).
     *   "id":       Unique identifier for this feature. Repeated identifiers throw an error. (Also accepts "name".) Required.
     *   "chrom":    Which chromosome it occurs on (Also accepts "chr" or "chromosome")
     *   "start":    Start position on the chromosome
     *   "stop":     Stop position on the chromosome (Also accepts "end")
     *   "position": Postion on chromosome (in place of "start" and "stop" for features that are a single nucleotide)
     *   "parent":   What other named feature this descends from (eg, Gene -> Transcript -> Exon). If none given, this
     *                   will default to the chromosome (or the genome if chromosome isn't supplied)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files).
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromJsonFile(String filename) {
        //TODO: Add JSON API and parse. Allow additional keys beyond just what's listed here in an extensible hash
        return this;
    }

    /**
     * Load in data from a flat, tab-delimited text file. The first row should be a header identifying what attribute is
     * in each column, and each subsequent row should correspond to a single feature. Columns that don't apply to a
     * given feature should use "NA" or left empty. Common attributes (columns) are listed below. Although only the "id"
     * attribute is required, a feature is pretty useless without some sort of positional information (chromosome,
     * start/stop, etc.).
     *   "id":       Unique identifier for this feature. Repeated identifiers throw an error (Also accepts "name".) Required.
     *   "chrom":    Which chromosome it occurs on (Also accepts "chr" or "chromosome")
     *   "start":    Start position on the chromosome
     *   "stop":     Stop position on the chromosome (Also accepts "end")
     *   "position": Postion on chromosome (in place of "start" and "stop" for features that are a single nucleotide)
     *   "parent":   What other named feature this descends from (eg, Gene -> Transcript -> Exon). If none given, this
     *                   will default to the chromosome (or the genome if chromosome isn't supplied)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files).
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromFlatFile(String filename) {
        //TODO: Load columns into hashes to make features extensible? - maybe for optional things
        return this;
    }

}
