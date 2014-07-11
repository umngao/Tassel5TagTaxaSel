package net.maizegenetics.dna.map;

import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import net.maizegenetics.util.DirectedGraph;
import net.maizegenetics.util.Utils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by jgw87 on 7/2/14.
 * A map to hold genome features for lookup by name and by location. The features themselves are hierarchical and so
 * can be traced up and down the tree.
 *
 * As methods are added to return based on different filters, try to unify the results. That is, all filter-like methods
 * should return the same sort of data structure, such as a HashSet of GenomeFeatures. It may be worthwhile to create
 * a custom GenomeFeatureSet class to be able to string such operations together (mygenes = mygenes.ofType("exon").onChrom(1).inRange(1, 10000);),
 * although whether such filters are needed for the bulk of this class's purpose (matching SNPs to genome annotations) is yet
 * to be seen.
 *
 * This class shouldn't be created directly, but should instead use the GenomeFeatureMapBuilder class .build() method
 */
//TODO: Add functionality to write out a full constructed map to a file of some sort
public class GenomeFeatureMap {

    //Root of the GenomeFeature tree - DEPRECATED. REPLACED BY AN EXPLICIT GRAPH, BELOW
    //GenomeFeature mygenome = null;

    DirectedGraph<GenomeFeature> featureTree = null;

    //Lookups to identify GenomeFeatures by their name
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private Multimap<String, GenomeFeature> typeLookup = null;
    //TODO: Figure out how to organize this, especially for dealing with multiple chromosomes. Nested lookup structure?
    private HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup = null;

    /**
     * Default constructor for creating a GenomeFeatureMap from pre-made data structures. This should ONLY be called from
     * the {@link GenomeFeatureMapBuilder} class
     * @param nameLookup Lookup table of unique IDs -> {@link GenomeFeature} objects
     * @param locationLookup Lookup table to retrieve {@link GenomeFeature}s by their genomic location
     * @param typeLookup Lookup table to retrieve {@link GenomeFeature}s by their type (exon, UTR, etc)
     * @param featureTree The graph of genomic features. Should be a directed graph, with the first two levels being genome and chromosome
     */
    GenomeFeatureMap(HashMap<String, GenomeFeature> nameLookup, Multimap<String, GenomeFeature> typeLookup,
                     HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup, DirectedGraph<GenomeFeature> featureTree){
        //this.mygenome=root;
        this.typeLookup=typeLookup;
        this.nameLookup=nameLookup;
        this.locationLookup=locationLookup;
        this.featureTree = featureTree;
    }

    public GenomeFeature getFeatureFromId(String id){
        return nameLookup.get(id);
    }

    public Collection<GenomeFeature> getFeaturesAtLocation(int chrom, int position){
        return getFeaturesInRange(chrom, position, position);
    }

    //TODO: Figure out how to get this working
    public Collection<GenomeFeature> getFeaturesInRange(int chrom, int start, int end){
        Range myrange = Range.closed(start, end); //'Closed' = inclusive, so closed(1,3) = 1,2,3 and closed(1,1) = 1
        return null;
    }


    /**
     * Write just the location lookup to a tab-delimited file. This is mostly to check that your locations loaded properly,
     * as there is no way to read them back in.
     * @param filename The output file to be written to
     */
    //TODO: Test and use this
    public void writeLocationLookupToFile(String filename){
        try {
            BufferedWriter writer = Utils.getBufferedWriter(filename);
            writer.append("chrom\tstart\tstop\tfeatures\n");
            String[] chroms =locationLookup.keySet().toArray(new String[0]);
            Arrays.sort(chroms);
            for(String chrom: chroms){
                Map<Range<Integer>, HashSet<GenomeFeature>> itermap = locationLookup.get(chrom).asMapOfRanges();
                for(Range<Integer> r: itermap.keySet()){
                    int start = r.lowerEndpoint();
                    int stop = r.upperEndpoint();
                    writer.append(chrom + "\t" + start + "\t" + stop + "\t");
                    //List of features
                    HashSet<GenomeFeature> features = itermap.get(r);
                    for(GenomeFeature f: features) {
                         writer.append(f.id() +";");
                    }
                    writer.append("\n");
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Write the map data as a JSON file (which can be read in by {@link GenomeFeatureMapBuilder}). Core attributes
     * (unique ID, chromosome, start, stop, and parent ID) are output to all features. Any additional attributes are
     * output only for those features that have them. Attributes are output in alphabetical order. Since the attribute
     * name has to be output for every feature, this can waste space if all your features have the same attributes. In
     * that case a tab-delimited flat file or precompiled binary file is probably a better choice.
     * @param filename The output file to be written to
     */
    //TODO: Test and use this
    public void writeMapAsJsonFile(String filename){

    }

    /**
     * Write the map data as a flat, tab-delimited file (which can be read in by {@link GenomeFeatureMapBuilder}). The
     * writer first compiles a set of all attribute types across the map and sets these as the columns (in alphabetical
     * order). Attributes that don't apply to a given feature are output as "NA". This can end up being very wasteful
     * if you have some attributes that only apply to a small percentage of features; in that case, a JSON or binary
     * file is probably a better choice.
     * @param filename The output file to be written to
     */
    //TODO: Test and use this
    public void writeMapAsFlatFile(String filename){

    }

    /**
     * Write the map data as a precompiled binary that can be read back in by {@link GenomeFeatureMapBuilder}). Unlike
     * JSON or flatfile format, this file is not human-readable and is simply a a representation of the Java data
     * structures saved onto a disk. This makes it very compact and fast to read back in. This is the preferred format
     * for long-term storage of a {@link GenomeFeatureMap} since there is (almost) no risk of someone accidentally modifying
     * the file.
     * @param filename The output file to be written to
     */
    //TODO: Test and use this
    //TODO: Implement Serializable interface in order to write all these things out
    public void writeMapAsBinaryFile(String filename){

    }
}
