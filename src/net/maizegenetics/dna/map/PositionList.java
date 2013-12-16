package net.maizegenetics.dna.map;

import net.maizegenetics.dna.snp.GenotypeTable;
import java.util.List;

/**
 * List of positions in the genome. This type is used by every
 * {@link GenotypeTable}, but it can also be used list of GWAS results and other
 * genomic annotations.
 *
 * @author Terry Casstevens and Ed Buckler
 */
public interface PositionList extends List<Position> {

    /**
     * Return reference diploid allele values at given site.
     *
     * @param site site
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte referenceGenotype(int site);

    /**
     * Returns reference sequence of diploid allele values for given taxon in
     * specified range (end site not included). Each value in array contains
     * both diploid values. First four bits holds the first allele, and the
     * second four bits holds the second allele.
     *
     * @param startSite start site
     * @param endSite end site
     *
     * @return reference sequence of diploid allele values.
     */
    public byte[] referenceGenotypes(int startSite, int endSite);

    /**
     * Returns reference sequence of diploid allele values. Each value in array
     * contains both diploid values. First four bits holds the first allele, and
     * the second four bits holds the second allele.
     *
     * @return reference sequence of diploid allele values.
     */
    public byte[] referenceGenotypeForAllSites();

    /**
     * Return whether this alignment has defined reference sequence.
     *
     * @return true if this alignment has reference sequence.
     */
    public boolean hasReference();

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String siteName(int site);

    /**
     * Returns total number of sites of this alignment.
     *
     * @return number of sites
     */
    public int numberOfSites();

    /**
     * Return number of sites for given Chromosome
     *
     * @param chromosome
     * @return number of sites
     */
    public int chromosomeSiteCount(Chromosome chromosome);

    /**
     * Get the first (inclusive) and last (exclusive) site of the specified
     * chromosome in this alignment.
     *
     * @param chromosome chromosome
     *
     * @return first and last site
     */
    public int[] startAndEndOfChromosome(Chromosome chromosome);

    /**
     * Returns the physical position at given site.
     *
     * @param site site
     *
     * @return physical position
     */
    public int chromosomalPosition(int site);

    /**
     * Return site of given physical position in chromosome. If the physical
     * position doesn't exist, (-(insertion point) - 1) is returned. If
     * chromosome is not found, an exception is thrown.
     *
     * @param physicalPosition physical position
     * @param chromosome chromosome. if null, the first chromosome is used.
     *
     * @return index
     */
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome);

    /**
     * Return site of given physical position / SNP ID in chromosome. If the
     * physical position doesn't exist, (-(insertion point) - 1) is returned. If
     * chromosome is not found, an exception is thrown. This is to support
     * multiple sites with the same physical position but different SNP IDs.
     *
     * @param physicalPosition physical position
     * @param chromosome chromosome. if null, the first chromosome is used.
     * @param snpName SNP ID
     *
     * @return index
     */
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName);

    /**
     * Returns all physical positions.
     *
     * @return physical positions.
     */
    public int[] physicalPositions();

    /**
     * Return Chromosome Name for given site.
     *
     * @param site site
     *
     * @return Chromosome Name
     */
    public String chromosomeName(int site);

    /**
     * Return Chromosome for given site.
     *
     * @param site site
     *
     * @return Chromosome
     */
    public Chromosome chromosome(int site);

    /**
     * Return Chromosome with matching name. First to match will be returned.
     *
     * @param name name
     *
     * @return Chromosome
     */
    public Chromosome chromosome(String name);

    /**
     * Return all chromosomes.
     *
     * @return chromosomes
     */
    public Chromosome[] chromosomes();

    /**
     * Return number of chromosomes.
     *
     * @return number of chromosomes
     */
    public int numChromosomes();

    /**
     * Returns starting site for each chromosome.
     *
     * @return starting site for each chromosome.
     */
    public int[] chromosomesOffsets();

    /**
     * Return size of indel at given site.
     *
     * @param site site
     *
     * @return indel size
     */
    public int indelSize(int site);

    /**
     * Returns whether give site is an indel.
     *
     * @param site site
     *
     * @return true if indel
     */
    public boolean isIndel(int site);

    /**
     * Gets the Genome Assembly.
     *
     * @return the genome assembly.
     */
    public String genomeVersion();

    /**
     * Return whether is positive strand at given site.
     *
     * @param site site
     *
     * @return whether is positive strand.
     */
    public boolean isPositiveStrand(int site);

}