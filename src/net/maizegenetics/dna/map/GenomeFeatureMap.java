package net.maizegenetics.dna.map;

import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

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

    //Root of the GenomeFeature tree
    GenomeFeature mygenome = null;

    //Lookups to identify GenomeFeatures by their name
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private Multimap<String, GenomeFeature> typeLookup = null;
    //TODO: Figure out how to organize this, especially for dealing with multiple chromosomes. Nested lookup structure?
    private RangeMap<Integer, HashSet<GenomeFeature>> locationLookup = null;

    /**
     * Default constructor for creating a GenomeFeatureMap from pre-made data structures. This should ONLY be called from
     * the GenomeFeatureMapBuilder class
     * @param root Root of the GenomeFeature tree
     * @param nameLookup Lookup table of unique IDs -> GenomeFeature objects
     * @param locationLookup Lookup table to retrieve GenomeFeatures by their genomic location
     */
    GenomeFeatureMap(GenomeFeature root, HashMap<String, GenomeFeature> nameLookup,
                     Multimap<String, GenomeFeature> typeLookup, RangeMap<Integer, HashSet<GenomeFeature>> locationLookup){
        this.mygenome=root;
        this.typeLookup=typeLookup;
        this.nameLookup=nameLookup;
        this.locationLookup=locationLookup;
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

}
