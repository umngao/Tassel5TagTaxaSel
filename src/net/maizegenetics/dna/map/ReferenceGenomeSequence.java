/**
 * 
 */
package net.maizegenetics.dna.map;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import net.maizegenetics.analysis.gbs.DiscoverySNPCallerPlugin;
import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Utils;

/**
 * ReferenceGenomeSequence class.  This class is used to read chromosome sequences
 * from fasta files.  Data is stored as half-bytes packed into a byte array.
 * This byte array comprises the "value" for a hash map whose key is a
 * Chromosome object.
 * 
 * The class also  methods to obtain a full or partial genome sequence for a
 * specified stored chromosome.
 * 
 * @author Lynn Johnson
 *
 */
public class ReferenceGenomeSequence implements GenomeSequence{

	private static final Logger myLogger = Logger.getLogger(ReferenceGenomeSequence.class);

	private Map<Chromosome, byte[]> chromPositionMap;

	ReferenceGenomeSequence(Map<Chromosome, byte[]>chromMap){
		this.chromPositionMap = chromMap;
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
		// This code written assuming caller sent in 1 based number

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
						return null;
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

                        if (!evenNumSites){ // add last byte
                            fullBytes[fIndex] = (byte) (packedBytes[pIndex] >> 4);
                            fullBytes[fIndex+1] = (byte)(packedBytes[pIndex] & 0x0F);
                        }
 
						if (evenNumSites) {
							// One more allele to store - top half of last byte 
	                       	fullBytes[fIndex] = (byte)(packedBytes[pIndex] >> 4); 
						}
						return fullBytes;
					} 
				}
			}
		}		
		return null;
	}

	protected static byte[] readReferenceGenomeChr(String fastaFileName, int targetChr) {
		// Read specified file, return entire sequence for requested chromosome 
		int nBases = 0;
		int basesPerByte = 2;

		myLogger.info("\n\nReading from the reference genome fasta file: " + fastaFileName);
		String line = null;
		try {
			boolean found = false;
			ArrayList<String> mySequence = new ArrayList<String>();
			BufferedReader br = Utils.getBufferedReader(fastaFileName); // this takes care of .gz
			int currChrom = Integer.MIN_VALUE;
			
			while ((line = br.readLine()) != null && !found) {
				line = line.trim();

				if (line.startsWith(">")) {
					String chrS = line.replace(">", "");
					chrS = chrS.replace("chromosome ", ""); // either chromosome, or chr or just number in file
					chrS = chrS.replace("chr", "");
					try {
						currChrom = Integer.parseInt(chrS);
					} catch (NumberFormatException e) {
						// If it fails, look for  keyFile data translation from DiscoverySNPCallerPlugin
						currChrom = DiscoverySNPCallerPluginV2.keyFileReturnChromInt(line);
						if (currChrom == -1) {
							myLogger.error("\n\nTagsToSNPByAlignment detected a non-numeric chromosome name in the reference genome sequence fasta file: " + chrS
									+ "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
									+ "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n"
									+ "Alternately, you may supply a tab-separated key file that maps each chromosome name to a number\n\n");
						}
					}
					myLogger.info("Currently reading chromosome " + currChrom +  ")");
				} else if (currChrom == targetChr) {
					mySequence.add(line);
					found = true;
					nBases += line.length();

					while ((line = br.readLine()) != null) {
						line = line.trim();

						if (line.startsWith(">")) {
							break; // end of chromosome sequence
						} else {
							// Store the bases for requested chromosome
							nBases += line.length();
							mySequence.add(line);
						}						
					}
				}
			}

			br.close();
			myLogger.info("\n\nFinished reading target chromosome " + targetChr + "\n");
			// Process the chromosome sequence data
			if (mySequence.size() > 0){
				int nBytes = (nBases % basesPerByte == 0) ? nBases/basesPerByte : (nBases / basesPerByte) + 1;
				byte[] sequence = new byte[nBytes];
				ByteBuffer bBuf = ByteBuffer.wrap(sequence);
				for (int sIndex = 0; sIndex < mySequence.size(); sIndex++) {
					String sLine = mySequence.get(sIndex);
					int currentStrPos = 0; 
					int byteBufferPos = 0;
					
					int nLineBytes = (sLine.length() % basesPerByte == 0) ? sLine.length()/basesPerByte : (sLine.length() / basesPerByte) + 1;
					byte[] lineSeq = new byte[nLineBytes];

					while (currentStrPos < sLine.length()-1) {
						// store encoded alleles as half bytes in the byte array
						byte halfByteUpper = NucleotideAlignmentConstants.getNucleotideAlleleByte(sLine.charAt(currentStrPos));
						byte halfByteLower = NucleotideAlignmentConstants.getNucleotideAlleleByte(sLine.charAt(currentStrPos+1));
						byte positionByte = (byte) ( (halfByteUpper << 4) | (halfByteLower));

						lineSeq[byteBufferPos] = positionByte;
						currentStrPos += 2;
						byteBufferPos++;
					}
										
					// Catch the last byte if line contained an odd number of bytes
					// if lines all have odd number of bytes, the end of each should be
					// stored with the beginning of the next one. Code assume file will
					// have even number of bytes except perhaps for last line of chromosome sequence
					if (currentStrPos < sLine.length()) {
						byte halfByteUpper = NucleotideAlignmentConstants.getNucleotideAlleleByte(sLine.charAt(currentStrPos));
						byte positionByte = (byte) (halfByteUpper << 4);
						lineSeq[byteBufferPos] = positionByte;
					}
					bBuf.put(lineSeq); // add to return buffer
				}
				return sequence;
			}

			br.close();
		} catch (Exception e) {
			myLogger.error("Exception caught while reading the reference genome fasta file at line.  Error=" + e);
			e.printStackTrace();
		}
		return null;
	}

}
