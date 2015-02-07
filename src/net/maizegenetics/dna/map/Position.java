package net.maizegenetics.dna.map;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.util.GeneralAnnotation;

/**
 * Defines a genomic positions and its known variants. Includes attributes of
 * chromosome, position, strand, centiMorgans, name (or SNP ID), whether this
 * position is a nucleotide, or includes an indel.
 *
 * @author Ed Buckler
 */
public interface Position extends Comparable<Position> {

    /**
     * Return the locus (generally a chromosome) of a site
     */
    public Chromosome getChromosome();

    /**
     * Return the physical position of a site
     */
    public int getPosition();

    /**
     * Return the strand for a site definition
     */
    public byte getStrand();

    /**
     * Return the strand for a site definition
     */
    public float getCM();

    /**
     * Return the ID (name) for a site
     */
    public String getSNPID();

    /**
     * Whether the position is a nucleotide position or another marker type
     * (SSR, AFLP, RAPD, CNV, which are recoded with text states)
     */
    public boolean isNucleotide();

    /**
     * Whether the position includes indels, which would be defined in the
     * variants
     */
    public boolean isIndel();

    /**
     * Returns the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or
     * {"100","103","106"}
     */
    public String[] getKnownVariants();

    /**
     * Return the minor allele frequency in a global scope
     */
    public float getGlobalMAF();

    /**
     * Returns the proportion of genotypes scored at a given site
     */
    public float getGlobalSiteCoverage();

    /**
     * Return the allele specified by alleleType, if unknown Alignment.Unknown
     * is return
     */
    public byte getAllele(WHICH_ALLELE alleleType);

    public GeneralAnnotation getAnnotation();

}
