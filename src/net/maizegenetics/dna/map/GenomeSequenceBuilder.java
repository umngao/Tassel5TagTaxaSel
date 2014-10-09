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
import java.util.stream.Collectors;

import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
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
	public GenomeSequence build() {
		return new HalfByteGenomeSequence(null);
	}

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
                        currChr=new Chromosome(currChr.getName(),currSeq.size(),null);
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
//                System.out.println(currSeq.size());
                currChr=new Chromosome(currChr.getName(),currSeq.size(),null);
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

		int fIndex = 0, pIndex = 0; // fullbyte index, packed byte index  		
		while (fIndex < unpkSequence.length -1) {
			byte halfByteUpper = NucleotideAlignmentConstants.getNucleotideAlleleByte((char)unpkSequence[fIndex]);
			byte halfByteLower = NucleotideAlignmentConstants.getNucleotideAlleleByte((char)unpkSequence[fIndex+1]);

			packedSequence[pIndex] = (byte) ( (halfByteUpper << 4) | (halfByteLower));
			fIndex +=2;
			pIndex++;
		}
		// Catch the last byte if unpkSequence contained an odd number of bytes
		if (fIndex < unpkSequence.length) {
			byte halfByteUpper = NucleotideAlignmentConstants.getNucleotideAlleleByte((char)unpkSequence[fIndex]);
			packedSequence[pIndex] = (byte) (halfByteUpper << 4);
		}
		return packedSequence;		
	}

	private static  Chromosome parseChromosome (String chromString) {
		int chromNumber = Integer.MAX_VALUE;
		String chrS = chromString.replace(">","");
		chrS = chrS.replace("chromosome ", ""); // either chromosome, or chr or just number in file
		chrS = chrS.replace("chr", "");
		Chromosome chromObject = null;
		try {
			chromNumber = Integer.parseInt(chrS);
			chromObject = new Chromosome(chrS);
		} catch (NumberFormatException e) {
			// If it fails, look for  keyFile data translation from DiscoverySNPCallerPlugin
			// This will be changed/deleted in favor of either (a) nothing or (b) data stored in database table
			chromNumber = DiscoverySNPCallerPluginV2.keyFileReturnChromInt(chromString);  
			if (chromNumber == -1) {
				myLogger.error("\n\nTagsToSNPByAlignment detected a non-numeric chromosome name in the reference genome sequence fasta file: " + chrS
						+ "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
						+ "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n");
			} else {
				chromObject = new Chromosome(Integer.toString(chromNumber));
			}
		}
		return chromObject;
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
        if (!chromPositionMap.isEmpty()) {
            return chromPositionMap.keySet();
        }
        return null;
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom) {
        return chromosomeSequence(chrom,1,chromLengthLookup.get(chrom));
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int lastSite) {
        startSite--;  //shift over to zero base
        lastSite--;   //shift over to zero base
        byte[] fullBytes = new byte[lastSite - startSite + 1];
        byte[] packedBytes = chromPositionMap.get(chrom);
        for (int i = startSite; i <= lastSite; i++) {
            fullBytes[i - startSite] = (byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F));
        }
        return fullBytes;
    }
}