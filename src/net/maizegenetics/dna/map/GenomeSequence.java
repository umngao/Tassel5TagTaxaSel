/**
 * Interface for Genome Sequence
 * 
 */
package net.maizegenetics.dna.map;

import java.util.Set;


/**
 * Defines the genome sequence of a chromosome
 * 
 * @author lcj34
 *
 */
public interface GenomeSequence {
	
	// Returns a set of chromosomes
	public Set<Chromosome> chromosomes();
    
    //Returns complete sequence for the specified chromosome    
    public byte[] chromosomeSequence(Chromosome chrom);
    
    //Returns a chromosome's genomic sequence contained by indicated start/end positions  
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int endSite);
    
    
    // Read fasta file in compressed (.gz) or uncompressed format
    public byte[] readReferenceGenomeChr (String fastaFileName, int targetChr);

}
