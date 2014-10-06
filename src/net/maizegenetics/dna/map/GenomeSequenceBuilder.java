/**
 * 
 */
package net.maizegenetics.dna.map;

import java.util.HashMap;
import java.util.Map;


/**
 * Builder for a chromosome genome sequence that hides the ReferenceGenomeSequence implementation.
 * The sequence is read from a fasta file and stored in a byte array, 2 alleles packed per byte.
 * 
 * THe intent is to move the ReferenceGenomeSequence class into this file.
 * 
 * @author Lynn Johnson
 *
 */
public class GenomeSequenceBuilder {
    private Map<Chromosome, byte[]> chromPositionMap = new HashMap<Chromosome, byte[]>();

    public GenomeSequenceBuilder( ) {
    }
 
    public GenomeSequence build() {
    	return new ReferenceGenomeSequence(chromPositionMap);
    }

    private GenomeSequenceBuilder(Chromosome myChrom, byte[] seqBytes) {
    	chromPositionMap.put(myChrom, seqBytes);
    }
    
    public static GenomeSequenceBuilder instance(String fastaFileName, int targetChr) {
    	byte[] seqBytes = ReferenceGenomeSequence.readReferenceGenomeChr(fastaFileName, targetChr);
       	if (seqBytes == null ) return null;
    	Chromosome myChrom = new Chromosome(String.valueOf(targetChr));
 
    	return new GenomeSequenceBuilder(myChrom,seqBytes);
    }
    
}