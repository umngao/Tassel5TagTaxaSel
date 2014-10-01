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
	Set<Chromosome> chromosomes();
    
    //Returns complete sequence for the specified chromosome    
    byte[] chromosomeSequence(Chromosome chrom);
    
    //Returns a chromosome's genomic sequence contained by indicated start/end positions  
    byte[] chromosomeSequence(Chromosome chrom, int startSite, int endSite);
    
    
    // Read fasta file in compressed (.gz) or uncompressed format
    byte[] readReferenceGenomeChr (String fastaFileName, int targetChr);

}
