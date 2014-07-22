package net.maizegenetics.dna.map;

import com.google.common.collect.*;
import net.maizegenetics.util.DirectedGraph;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import sun.org.mozilla.javascript.internal.json.JsonParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by jgw87 on 7/2/14.
 * A Builder class to create a GenomeFeatureMap to identify genomic features. Can build piecemeal or read in from a file
 * For now, not implementing functionality to make a new builder from an existing map.
 */
//TODO: Add functionality for reading in a precompiled map
//TODO: Add ability to read in GFF, Flatfile, and JSON
public class GenomeFeatureMapBuilder {

    //Root of the GenomeFeature tree - DEPRECATED IN FAVOR OF AN ACTUAL GRAPH
    //GenomeFeature root = null;

    private static final Logger myLogger = Logger.getLogger(GenomeFeatureMapBuilder.class);

    //Graph of all the genome features connecting to each other
    DirectedGraph<GenomeFeature> featureTree = null;

    //Lookups to identify GenomeFeatures by their name, location, and type (exon, gene, etc)
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup = new HashMap<>(); // Complex.
    private Multimap<String, GenomeFeature> typeLookup = HashMultimap.create();

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


        }

        //TODO: Remove 1-bp ranges caused by adjacent features; need a clean() function to call at end?
        //TODO: Or maybe add new features in a different way, finding and filling gaps rather than overlaying
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

        //Make ranges closed-open b/c otherwise results in a lot of 1-bp redundant ranges where adjacent features were added
        Range<Integer> featureRange = Range.closedOpen(feature.start(), feature.stop() + 1);

        //First, extract out subrange and save any annotations already there. (Requires making a copy so later modification of
        // masterMap doesn't change things
        ArrayList<Range<Integer>> rangeList = new ArrayList<>();
        ArrayList<HashSet<GenomeFeature>> hashList = new ArrayList<>();
        Map<Range<Integer>, HashSet<GenomeFeature>> subranges = masterMap.subRangeMap(featureRange).asMapOfRanges();
        for(Range<Integer> r: subranges.keySet()){
            rangeList.add(Range.closedOpen(r.lowerEndpoint(), r.upperEndpoint()));
            hashList.add(new HashSet<GenomeFeature>(subranges.get(r)));
        }



        //Next, set entire range equal to the new feature to cover any areas not covered by the existing ranges
        HashSet<GenomeFeature> newset = new HashSet<>();
        newset.add(feature);
        masterMap.put(featureRange, newset);

        //Take the extracted, existing annotations and add back in, adding this feature to the list since they overwrite
        // everything already there
        for(int i=0; i<rangeList.size(); i++){
            HashSet<GenomeFeature> tempset = hashList.get(i);
            Range<Integer> temprange = rangeList.get(i);
            tempset.add(feature); //Add the feature on top of existing ones
            masterMap.put(temprange, tempset);
        }
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
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files)
     *
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromGffFile(String filename){
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
     * (For example: "chromosome":1). Note that if you have more than one feature per file (the normal case), all but the
     * last closing brace ('}') should be followed by a comma, and the whole group should be within square braces ('[...]'
     * That is, the first character of the file should be '[' and the last should be ']'). This makes it a properly-formatted
     * JSON array.
     *
     * Common attributes (keys) are listed below. Although only the "id" attribute is required, a feature is pretty
     * useless without some sort of positional information (chromosome, start/stop, etc.).
     *   "id":       Unique identifier for this feature. Repeated identifiers throw an error. (Also accepts "name".) Required.
     *   "chrom":    Which chromosome it occurs on (Also accepts "chr" or "chromosome")
     *   "start":    Start position on the chromosome
     *   "stop":     Stop position on the chromosome (Also accepts "end")
     *   "position": Postion on chromosome (in place of "start" and "stop" for features that are a single nucleotide)
     *   "parent_id":   What other named feature this descends from (eg, Gene -> Transcript -> Exon). If none given, this
     *                   will default to the chromosome (or the genome if chromosome isn't supplied)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files).
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromJsonFile(String filename) {
        JSONParser jparse = new JSONParser();
        BufferedReader reader = Utils.getBufferedReader(filename);
        try {
            JSONArray jarray = (JSONArray) jparse.parse(reader);
            Iterator iter = jarray.iterator();
            while(iter.hasNext()){
                JSONObject json = (JSONObject) iter.next();
                GenomeFeature newFeature = new GenomeFeatureBuilder().parseJsonObject(json).build();
                addFeature(newFeature);
            }
            reader.close();
        } catch (IOException e) {
            myLogger.error("Error loading data from JSON file " + filename);
            e.printStackTrace();
        } catch (ParseException e) {
            myLogger.error("Error parsing information in JSON file " + filename);
            e.printStackTrace();
        }

        return this;
    }

    //Alternate implementation with manually reading in the JSON data; not recommended
    /*public GenomeFeatureMapBuilder addFromJsonFile(String filename) {
        JSONParser jparse = new JSONParser();
        BufferedReader reader = Utils.getBufferedReader(filename);
        try {
            String inline = reader.readLine();
            String tempJson = "";
            while(inline != null){
                tempJson += inline;
                //If has closing brace (mark of end of object), then parse
                while(tempJson.contains("}")){  //While loop in case multiple objects on one line
                    //Subtract out JSON object
                    int start=tempJson.indexOf('{');
                    int end=tempJson.indexOf('}');
                    String json = tempJson.substring(start, end+1);
                    System.out.println("jsonAsString:" + tempJson);
                    System.out.println("\ttake out:" + json);
                    tempJson = tempJson.substring(end + 1, tempJson.length());    //Remove everything
                    System.out.println("\tleft:" + tempJson);

                    //Add data as feature
                    JSONObject featureData = (JSONObject) jparse.parse(json);
                    GenomeFeature newFeature = new GenomeFeatureBuilder().parseJsonObject(featureData).build();
                    addFeature(newFeature);
                }
                inline = reader.readLine();
            }
        } catch (IOException e) {
            myLogger.error("Error loading data from JSON file " + filename);
            e.printStackTrace();
        } catch (ParseException e) {
            myLogger.error("Error parsing information in JSON file " + filename);
            e.printStackTrace();
        }

        return this;
    }*/

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
     *   "parent_id":   What other named feature this descends from (eg, Gene -> Transcript -> Exon). If none given, this
     *                   will default to the chromosome (or the genome if chromosome isn't supplied)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files).
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromFlatFile(String filename) {
        try {
            BufferedReader reader = Utils.getBufferedReader(filename);
            String line=reader.readLine();
            String[] header = null;
            int n=1;    //Keep track of line numbers
            while(line != null ){
                n++;
                if(line.startsWith("#")){    //Skip comment lines
                    line=reader.readLine();
                    continue;
                }
                String[] tokens = line.split("\t");

                if(header == null){ //Save header data
                    header=tokens;
                    line=reader.readLine();
                    continue;
                }

                //Check that number of fields is correct
                if(tokens.length != header.length){
                    myLogger.error("Error: line " + n + " has a different number of fields (" + tokens.length + ") than the header (" + header.length + ")");
                }

                //Load everything into a hapmap
                HashMap<String, String> data = new HashMap<>();
                for(int i=0; i<tokens.length; i++){
                    data.put(header[i], tokens[i]);
                }

                //Add to map
                GenomeFeature newFeature = new GenomeFeatureBuilder().loadAll(data).build();
                addFeature(newFeature);
                line=reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return this;
    }

}
