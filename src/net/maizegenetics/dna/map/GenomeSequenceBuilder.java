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
	private static Map<Chromosome, byte[]> chromPositionMap;

	public GenomeSequence build() {
		return new ReferenceGenomeSequence();
	}

	private GenomeSequenceBuilder(Map<Chromosome, byte[]> cMap) { 
		GenomeSequenceBuilder.chromPositionMap = cMap;
	}

	public static GenomeSequenceBuilder instance(String fastaFileName) {
		Map<Chromosome, byte[]> chromPositionMap = ReferenceGenomeSequence.readReferenceGenomeChr(fastaFileName);   	
		return new GenomeSequenceBuilder(chromPositionMap);
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
	protected static class ReferenceGenomeSequence implements GenomeSequence{
		private static final Logger myLogger = Logger.getLogger(ReferenceGenomeSequence.class);

		@Override
		public Set<Chromosome> chromosomes() {
			if (!chromPositionMap.isEmpty()) {
				return chromPositionMap.keySet();
			}
			return null;
		}

		@Override
		public byte[] chromosomeSequence(Chromosome chrom) {
			// Return all genome sequence bytes for specified chromosome
			int chromNumber = chrom.getChromosomeNumber();

			for (Map.Entry<Chromosome, byte[]> chromEntry : chromPositionMap.entrySet()) {
				int mapChromeNo = chromEntry.getKey().getChromosomeNumber();
				if (mapChromeNo == chromNumber) {
					byte[] packedBytes = chromEntry.getValue();
					int fullBytesLength = packedBytes.length * 2;
					byte[] fullBytes = new byte[fullBytesLength];
					int fIndex = 0, pIndex = 0;

					while (fIndex < fullBytesLength - 1){
						fullBytes[fIndex] = (byte) (packedBytes[pIndex] >> 4);
						fullBytes[fIndex+1] = (byte)(packedBytes[pIndex] & 0x0F);
						fIndex += 2;
						pIndex++;					
					}	
					return fullBytes;
				}
			}
			return null;
		}

		@Override
		public byte[] chromosomeSequence(Chromosome chrom, int startSite, int endSite) {
			// Return specified portion of a chromosome's sequence	
			// Code assumes caller sent in 1 based number

			if (startSite < 1 || startSite > endSite)
				return null;

			int chromNumber = chrom.getChromosomeNumber();
			int numSites = endSite - startSite + 1;  // start/end site values are inclusive
			boolean evenStart = startSite % 2 == 0 ? true: false;
			boolean evenNumSites = numSites % 2 == 0 ? true: false;

			byte[] fullBytes = null;
			if (numSites > 0) {			
				fullBytes = new byte[numSites];

				for (Map.Entry<Chromosome, byte[]> chromEntry : chromPositionMap.entrySet()) {
					int mapChromNo = chromEntry.getKey().getChromosomeNumber();
					if (mapChromNo == chromNumber) {
						byte[] packedBytes = chromEntry.getValue();

						if (startSite > packedBytes.length*2 ||
								endSite > packedBytes.length*2 ) {
							return null; // requested sequence is out of range
						}

						if (!evenStart) {
							// Start site is odd, we're grabbing a full byte to start
							// ex: user wants to start at position 5 in sequence.
							// sequence is stored as halfbyte, position 5 is the top
							// half of byte 2 in a 0-based byte array.  So offset is 5/2 = 2
							int fIndex = 0; // index into new full bytes array
							int pIndex = startSite/2; // index into stored pack bytes array

							while (fIndex < numSites - 1){
								fullBytes[fIndex] = (byte) (packedBytes[pIndex] >> 4);
								fullBytes[fIndex+1] = (byte)(packedBytes[pIndex] & 0x0F);
								fIndex += 2;
								pIndex++;                               
							}         
							if (numSites %2 != 0) {
								fullBytes[fIndex] = (byte)((packedBytes[pIndex] & 0xF0) >> 4);                        	
							}
							return fullBytes;
						} else {
							// start site is even number, we're starting in the lower half of a stored byte
							// ex: user wants to pull from position 4 in the sequence
							// Position 4 is stored in the lower half of byte 1 in the 0-based byte array
							// 4/2 - 1 = 1
							int pIndex = startSite/2 -1;
							int fIndex = 0;

							fullBytes[fIndex] = (byte)(packedBytes[pIndex] & 0x0F);
							fIndex++; pIndex++;
							while (fIndex < numSites - 2){
								fullBytes[fIndex] = (byte) (packedBytes[pIndex] >> 4);
								fullBytes[fIndex+1] = (byte)(packedBytes[pIndex] & 0x0F);
								fIndex += 2;
								pIndex++;                               
							}       
							if (evenNumSites) {
								fullBytes[fIndex] = (byte)(packedBytes[pIndex] >> 4); // store last allele, top half last byte
							} else { // store last byte
								fullBytes[fIndex] = (byte) (packedBytes[pIndex] >> 4); 
								fullBytes[fIndex+1] = (byte)(packedBytes[pIndex] & 0x0F);
							}
							return fullBytes;
						} 
					}
				}
			}		
			return null;
		}

		protected static  Map<Chromosome, byte[]> readReferenceGenomeChr(String fastaFileName) {
			// Read specified file, return entire sequence for requested chromosome 
			Map<Chromosome, byte[]> chromPositionMap = new HashMap<Chromosome, byte[]>();
			myLogger.info("\n\nReading from the reference genome fasta file: " + fastaFileName);

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
							chromPositionMap.put(currChr, halfByteCompression(currSeq.toByteArray()));
							myLogger.info("Finished reading chromosome " + currChr.getChromosomeNumber() +  ")");
						}
						currChr = parseChromosome(line); 
						currSeq = new ByteArrayOutputStream();
					} else {
						currSeq.write(line.getBytes());
					}
				}
				// reached end of file - write last bytes
				if (currSeq.size() > 0) {
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

			int nBytes = (unpkSequence.length % 2 == 0) ? unpkSequence.length/2 : (unpkSequence.length/2) + 1;
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
}