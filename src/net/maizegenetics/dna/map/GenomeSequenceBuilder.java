/**
 * 
 */
package net.maizegenetics.dna.map;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;


/**
 * Builder for a chromosome genome sequence that hides the ReferenceGenomeSequence implementation.
 * The sequence is read from a fasta file and stored in a byte array, 2 alleles packed per byte.
 * 
 * @author Lynn Johnson
 *
 */
public class GenomeSequenceBuilder {
	private static final Logger myLogger = Logger.getLogger(GenomeSequenceBuilder.class);
	
	public static GenomeSequence instance(String fastaFileName) {
		Map<Chromosome, byte[]> chromPositionMap = readReferenceGenomeChr(fastaFileName);   	
		return new HalfByteGenomeSequence(chromPositionMap);
	}

	protected static  Map<Chromosome, byte[]> readReferenceGenomeChr(String fastaFileName) {
		// Read specified file, return entire sequence for requested chromosome 
		Map<Chromosome, byte[]> chromPositionMap = new HashMap<Chromosome, byte[]>();
		Chromosome currChr = null;	
		ByteArrayOutputStream currSeq = new ByteArrayOutputStream();
		String line = null;
		try {
			boolean found = false;
			BufferedReader br = Utils.getBufferedReader(fastaFileName); // this takes care of .gz

			while ((line = br.readLine()) != null && !found) {
				line = line.trim();

				if (line.startsWith(">")) {
					if (currChr != null) {
						// end processing current chromosome sequence
                        currChr=new Chromosome(currChr.getName(),currSeq.size(),currChr.getAnnotation());
						chromPositionMap.put(currChr, halfByteCompression(currSeq.toByteArray()));
					}
					currChr = parseChromosome(line); 
					currSeq = new ByteArrayOutputStream();
				} else {
					currSeq.write(line.getBytes());
				}
			}
			// reached end of file - write last bytes
			if (currSeq.size() > 0) {
                currChr=new Chromosome(currChr.getName(),currSeq.size(),currChr.getAnnotation());
				chromPositionMap.put(currChr, halfByteCompression(currSeq.toByteArray()));
			}
			br.close();
		} catch (IOException ioe) {
			System.out.println("ReferenceGenomeSequence: caught buffered read exception");
		}
		return chromPositionMap;
	}

    private static byte[] halfByteCompression(byte[] unpkSequence){
        // Take byte array, turn bytes into NucleotideAlignmentConstant
        // allele values, store as half bytes
        int nBytes = (unpkSequence.length+1)/2;
        byte[] packedSequence = new byte[nBytes];
        for (int i = 0; i < unpkSequence.length; i++) {
            byte halfByte = NucleotideAlignmentConstants.getNucleotideAlleleByte((char)unpkSequence[i]);
            if(i%2==0) halfByte<<=4;
            packedSequence[i/2]|=halfByte;
        }
        return packedSequence;
    }

	private static  Chromosome parseChromosome (String chromString) {
		String chrS = chromString.replace(">","");
		chrS = chrS.replace("chromosome ", ""); // either chromosome, or chr or just number in file
		chrS = chrS.replace("chr", "");
		GeneralAnnotation myAnnotations = null;
		
		String currChrDesc = null;
		int spaceIndex = chrS.indexOf(" ");
		if (spaceIndex > 0) {			
			currChrDesc = chrS.substring(chrS.indexOf(" ") + 1);
			myAnnotations = GeneralAnnotationStorage.getBuilder().addAnnotation("Description", currChrDesc).build();
			chrS = chrS.substring(0,chrS.indexOf(" "));
		} 
		return new Chromosome(chrS, -1, myAnnotations);
	}
   
}

/**
 * ReferenceGenomeSequence class.  This class is used to read chromosome sequences
 * from fasta files.  Data is stored as half-bytes packed into a byte array.
 * This byte array comprises the "value" for a hash map whose key is a
 * Chromosome object.
 *
 * The class also contains methods to obtain a full or partial genome sequence for a
 * specified stored chromosome.
 *
 * @author Lynn Johnson
 *
 */
class HalfByteGenomeSequence implements GenomeSequence{
    private Map<Chromosome, byte[]> chromPositionMap;
    private Map<Chromosome, Integer> chromLengthLookup=new HashMap<>();


    protected HalfByteGenomeSequence(Map<Chromosome, byte[]>chromPositionMap) {
        this.chromPositionMap = chromPositionMap;
        chromPositionMap.entrySet().stream()
                .forEach(e -> chromLengthLookup.put(e.getKey(),e.getKey().getLength()));

    }
    @Override
    public Set<Chromosome> chromosomes() {
        return chromPositionMap.keySet();
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom) {
        return chromosomeSequence(chrom,1,chromLengthLookup.get(chrom));
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int lastSite) {
        startSite--;  //shift over to zero base
        lastSite--;   //shift over to zero base
        byte[] packedBytes = chromPositionMap.get(chrom);
        if (startSite > packedBytes.length*2 || lastSite > packedBytes.length*2 ) {
        	return null; // requested sequence is out of range
        }
        byte[] fullBytes = new byte[lastSite - startSite + 1];
        for (int i = startSite; i <= lastSite; i++) {
            fullBytes[i - startSite] = (byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F));
        }
        return fullBytes;
    }
}