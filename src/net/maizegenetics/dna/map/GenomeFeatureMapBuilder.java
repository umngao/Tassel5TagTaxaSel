package net.maizegenetics.dna.map;

import com.google.common.collect.Multimap;
import com.google.common.collect.RangeMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by jgw87 on 7/2/14.
 * A Builder class to create a GenomeFeatureMap to identify genomic features. Can build piecemeal or read in from a file
 * For now, not implementing functionality to make a new builder from an existing map.
 */
//TODO: Add functionality for reading in a precompiled map
public class GenomeFeatureMapBuilder {

    //Root of the GenomeFeature tree
    GenomeFeature root = null;

    //Lookups to identify GenomeFeatures by their name, location, and type (exon, gene, etc)
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private RangeMap<Integer, HashSet<GenomeFeature>> locationLookup = null;
    private Multimap<String, GenomeFeature> typeLookup = null;


    public GenomeFeatureMap build(){
        buildGenomeTree();
        buildLocationLookup();
        return new GenomeFeatureMap(root, nameLookup, typeLookup, locationLookup);
    }


    private void buildGenomeTree(){
        //TODO: Implement

        //TODO: Build chromosomes first

        //Pull out the "genome" feature to serve as root. If it doesn't exist, create it.
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
        }
    }

    private void buildLocationLookup(){
        //TODO: Implement
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
        return this;
    }


    //TODO: Add method to read in a GFF or (better) GTF file and populate everything automatically


}
