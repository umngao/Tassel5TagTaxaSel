package net.maizegenetics.dna.tag;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;

import java.util.Map;
import java.util.Set;

/**
 * Interface to GBS Tag Data.  Replaces TagsCounts, TBT, and TOPM
 * TODO: TAS-476
 * @author Ed Buckler
 */
public interface TagData {


    /**
     * Returns the TaxaDistribution for a given taxa.  If distribution data is not available, null is returned.
     */
    TaxaDistribution getTaxaDistribution(Tag tag);

    /**
     * Provides all SNP allele calls associated with a given tag.
     * @param tag used for query
     * @return Map of position with allele state.  The positions are the positions of the actual SNP
     * (not the start of the tag).  The Byte value is defined by {@link net.maizegenetics.dna.snp.NucleotideAlignmentConstants}.
     * If no allele are available for a tag, the map will be empty.
     */
    Set<Allele> getAlleles(Tag tag);


    /**
     * Provides an iterator that goes through all tags with non-empty alleleMaps.
     * @return an iterator of alleleMaps (Map of position with allele state.  The positions are the positions of the actual SNP
     * (not the start of the tag).  The Byte value is defined by {@link net.maizegenetics.dna.snp.NucleotideAlignmentConstants}.)
     */
    Multimap<Tag,Allele> getAlleleMap();


    /**
     * Return all tags that cover a particular SNP position.  Empty list is returned if there are no tags.
     */
    Set<Tag> getTagsForSNPPosition(Position position);

    /**
     * Return all tags that call a particular SNP position allele.  Empty list is returned if there are no tags.
     */
    Set<Tag> getTagsForAllele(Position position, byte allele);

    /**
     * Return all tags that call a particular SNP position allele.  Empty list is returned if there are no tags.
     */
    Set<Tag> getTagsForAllele(Allele allele);

    /**
     * Iterator to go through all tags
     */
    Set<Tag> getTags();

    /**
     * Create the unique list of SNPs (variants) based on all alleles.
     * @return PositionList of all SNPs
     */
    PositionList getSNPPositions();

    /**
     * Create the unique list of SNPs (variants) based on all alleles.
     * Only alleles with support values greater than
     * the minSupportValue are considered.
     * @return PositionList of all SNPs
     */
    PositionList getSNPPositions(int minSupportValue);
    
    /**
     * Get SNPs for specified chromosomes
     * @param starting chromosome number
     * @param ending chromosome number
     * @return  A multimap with Chromosome Number (integer) as the key and a list of SNP positions as the value
     */
    ListMultimap<Integer,Integer> getSNPPositionsForChromosomes(Integer startChr,Integer endChr);

    /**
     * Returns the list of tags present of a taxon.  Note this could be a very compute intensive request.
     */
    Set<Tag> getTagsForTaxon (Taxon taxon);


    /*
     * Return a map of Tags (keys) and depth (value) for a given taxon and SNP position.
     */
    Map<Tag, Integer> getTagDepth(Taxon taxon, Position position);

    /*
     * Return a map of Tags (keys) and depth (value) greater or equal to a minimum depth
     */
    Map<Tag, Integer> getTagsWithDepth(int minimumDepth);


    /**
     * Get all of the genomic cut positions associated with all tags.  Positions in the position list are annotated with
     * cigarAlignment, isBest, alignmentApproach, and supportValue.
     * @return iterator of annotated positions
     */
    PositionList getTagCutPositions(boolean onlyBest);

    /**
     * Get the unique list of position that tags to for a specific chromosome within a range.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @param onlyBest only return positions flagged as best
     * @return List of positions
     */
    PositionList getTagCutPositions(Chromosome chromosome, int firstPosition, int lastPosition, boolean onlyBest);

    /*
    Map of alignment approaches names (key) and their protocols (value)
     */
    Map<String,String> getTagAlignmentApproaches();

    /**
     * Map of positions and with associated map of Tags and their taxa distribution and their alignment direction.
     * Warning:  This can be a very large data structure for entire chromosomes. Only the best positions are returned.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @return  Map of maps for position (key) to Map of Tag(key) to the Tuple(Direction,TaxaDistribution)(Value)
     */
    Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> getCutPositionTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition);


    /**
     * For a given genomic position returns of the map tags and their distribution.
     * @param cutPosition Genomic position of the cut site.
     * @return Map of Tag(key) TaxaDistribution(Value)
     */
    Map<Tag,TaxaDistribution> getTagsTaxaMap(Position cutPosition);




}
