/**
 * Interface for Genome Sequence
 * 
 */
package net.maizegenetics.dna.map;

import java.util.Set;


/**
 * Defines the genome sequence of a chromosome
 * 
 * @author Lynn Johnson
 *
 */
public interface GenomeSequence {
	
	/**
	 * Returns a list of chromosomes whose sequences have been
	 * stored in the chromsomeSequence map of the class implementing
	 * this interface.
	 * 
	 * @return  a Set of Chromosome objects
	 */
	public Set<Chromosome> chromosomes();
    
	/**
	 * Takes a Chromosome object and returns the stored byte array representing
	 * the genomic sequence for the specified chromosome.
	 * 
	 * @param chrom; a Chromosome object representing the chromosome whose
	 * 				sequence will be returned
	 * @return  A byte array containing the chromosome alleles in NucleotideAlignmentConstant
	 * 			form packed in half bytes
	 */
    public byte[] chromosomeSequence(Chromosome chrom);
    
    /**
     * Returns the partial genomic sequence for a  chromosome, from the specified start
     * position to the specified end position.  THe start/end positions are inclusive and
     * the request is 1-based (though the alleles are stored in a 0-based byte array).
     * 
     * @param chrom:  the chromosome whose partial sequence will be returned.
     * @param startSite:  the 1-based position in the sequence to start the pull.
     * @param endSite:  the 1-based position in the sequence that will be the last allele in the pull
     * @return A byte array of alleles in NucleotideAlignmentConstant form that is packed into
     * 			half bytes.
     */
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int endSite);

}
