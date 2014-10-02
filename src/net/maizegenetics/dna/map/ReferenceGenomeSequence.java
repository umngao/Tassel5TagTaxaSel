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
 * This byte array comprises the "value" for a hash map whose key is an int
 * indicating a chromosome.
 * 
 * The class also provides a Java Set that contains a list of chromosomes
 * for which we have stored a sequence, and methods to obtain a full or
 * partial genome sequence for a specified stored chromosome.
 * 
 * @author Lynn Johnson
 *
 */
public class ReferenceGenomeSequence implements GenomeSequence{

	private static final Logger myLogger = Logger.getLogger(ReferenceGenomeSequence.class);

	// Map to store chromosomes and sequence bytes 	
	private Map<String, byte[]> chromPositionMap = new LinkedHashMap<String, byte[]>();

	// Using HashSet for performance.  If we need an order maintained, this should be
	// changed to TreeSet
	private Set<Chromosome> chromosomeList = new HashSet<Chromosome>();

	@Override
	public Set<Chromosome> chromosomes() {
		return chromosomeList;
	}

	@Override
	public byte[] chromosomeSequence(Chromosome chrom) {
		// Return all genome sequence bytes for specified chromosome
		int chromNumber = chrom.getChromosomeNumber();

		for (Map.Entry<String, byte[]> chromEntry : chromPositionMap.entrySet()) {
			int mapChromeNo = Integer.parseInt(chromEntry.getKey());
			if (mapChromeNo == chromNumber) {
				return chromEntry.getValue();
			}
		}
		return null;
	}

	@Override
	public byte[] chromosomeSequence(Chromosome chrom, int startSite, int endSite) {
		// Return specified portion of a chromosome's sequence	
		// This code written assuming caller sent in 1 based number

		int chromNumber = chrom.getChromosomeNumber();
		int numSites = endSite - startSite + 1;  // assume start/end site values are inclusive
		boolean evenStart = startSite % 2 == 0 ? true: false;
		boolean evenNumSites = numSites % 2 == 0 ? true: false;
		int nBytes = numSites/2;
		if (evenStart) {
			nBytes++;
		} else if (!evenNumSites)  {
			nBytes++;			
		}

		byte[] sequence = null;
		if (nBytes > 0) {			
			sequence = new byte[nBytes];
			ByteBuffer bBuf = ByteBuffer.wrap(sequence);

			for (Map.Entry<String, byte[]> chromEntry : chromPositionMap.entrySet()) {
				int mapChromNo = Integer.parseInt(chromEntry.getKey());

				if (mapChromNo == chromNumber) {
					byte[] chromBytes = chromEntry.getValue();
					int seqLength = chromBytes.length;

					if (!evenStart) {
						// Start site is odd, we'll grabbing a full byte to start
						// ex: user wants to start at position 5 in sequence.
						// sequence is stored as halfbyte, position 5 is the top
						// half of byte 2 in a 0-based byte array.  So offset is 5/2 = 2
						int offset = startSite / 2;  // For 0 based array, this gets us to correct beginning entry
						int endPull = offset + numSites;

						// Verify range is valid
						if (( offset < seqLength) && ( endPull <= seqLength ))  {
							bBuf.put(chromEntry.getValue(),offset,nBytes-1);
							// last byte is either half byte or whole, depending on an even or odd number of bytes
							if (numSites %2 == 0) {
								// add the full last byte
								bBuf.put(chromBytes, offset+nBytes-1,1);
							} else {
								// we only want the top half of the last byte 
								byte lastByte = chromBytes[offset+nBytes-1];
								byte byteToStore = (byte) ((lastByte & 0xF0));
								bBuf.put(byteToStore);
							}
						} else {
							System.out.println("ReferenceGenomeSequence:chromosomeSequence: requested bytes are out of range for Chromosome");
							return null;
						}
						return sequence;

					} else {
						// start site is even number, we're starting in the lower half of a stored byte
						// ex: user wants to pull from position 4 in the sequence
						// Position 4 is stored in the lower half of byte 1 in the 0-based byte array
						// 4/2 - 1 = 1
						// All bytes have to be shifted for the return array
						int offset = startSite / 2 - 1;
						int endPull = offset + numSites;
						if (( offset < seqLength) && ( endPull <= seqLength ))  {
							byte[] firstByteArray = new byte[nBytes];
							ByteBuffer b2Buf = ByteBuffer.wrap(firstByteArray);
							b2Buf.put(chromBytes, offset,nBytes);

							int shiftedNBytes = evenNumSites == true ? nBytes-1 : nBytes;

							byte[] shiftedByteArray = new byte[shiftedNBytes]; // source array has one extra byte
							ByteBuffer b3Buf = ByteBuffer.wrap(shiftedByteArray);

							// get first byte, then loop
							int index= 0;

							byte firstByte = firstByteArray[index];
							byte secondByte = 0;
							byte tempByte = 0;
							for (index=1; index < nBytes; index++) {
								// Grab bytes and shift.  Lower half of first byte and upper half
								// of second byte make up a byte for the shifted array
								secondByte = firstByteArray[index];
								byte shiftArrayByte;
								if (index == 1) {
									shiftArrayByte = (byte) (((firstByte & 0x0F) << 4) | ((secondByte & 0xF0) >> 4));
								} else {
									shiftArrayByte = (byte) ((tempByte << 4) | ((secondByte & 0xF0) >> 4));
								}

								b3Buf.put(shiftArrayByte); 
								tempByte = (byte) (secondByte & 0x0F);	// process lower 4 bits of the byte									
							}
							if (!evenNumSites) {
								// One more byte to store
								byte shiftArrayByte = (byte)(secondByte << 4);
								b3Buf.put(shiftArrayByte);
							}
							return shiftedByteArray;
						} else {
							System.out.println("ReferenceGenomeSequence:chromosomeSequence: requested bytes are out of range for Chromosome");
						}
					} // end offset validation					
				}
			}
		}		
		return null;
	}

	@Override
	public  byte[] readReferenceGenomeChr(String fastaFileName, int targetChr) {
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
						byte halfByteUpper = getAlleleEncoding(sLine.charAt(currentStrPos));						
						byte halfByteLower = getAlleleEncoding(sLine.charAt(currentStrPos+1));
						byte positionByte = (byte) ( (halfByteUpper << 4) | (halfByteLower));

						lineSeq[byteBufferPos] = positionByte;
						currentStrPos += 2;
						byteBufferPos++;
					}
										
					// Catch the last byte if line contained an odd number of bytes
					// if lines all have odd number of bytes, the end of each should be
					// stored with the beginning of the next one.  I'm guessing we will
					// have even number of bytes except perhaps for last line of chromosome sequence
					if (currentStrPos < sLine.length()) {
						byte halfByteUpper = getAlleleEncoding(sLine.charAt(currentStrPos));
						byte positionByte = (byte) (halfByteUpper << 4);
						lineSeq[byteBufferPos] = positionByte;
					}
					bBuf.put(lineSeq); // add to return buffer
				}

				// Add chromosome to Set of Chromosomes.  
				Chromosome chrObj = new Chromosome(String.valueOf(currChrom));
				if (!chromosomeList.contains(chrObj)) {
					chromosomeList.add(chrObj);
				}
				String chromNumber = Integer.toString(currChrom);
				chromPositionMap.put(chromNumber, sequence); // store chromosome and bytes in the map			
				return sequence;
			}

			br.close();
		} catch (Exception e) {
			myLogger.error("Exception caught while reading the reference genome fasta file at line.  Error=" + e);
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * returns NucleotideAlignmentConstant encoding for a particular byte
	 */
	private byte getAlleleEncoding(char myChar) {
		byte myByte = 0x00;

		switch (myChar) {
		case 'A':
		case 'a':
			myByte = NucleotideAlignmentConstants.A_ALLELE;
			break;
		case 'C':
		case 'c':
			myByte = NucleotideAlignmentConstants.C_ALLELE;
			break;
		case 'G':
		case 'g':
			myByte = NucleotideAlignmentConstants.G_ALLELE;
			break;
		case 'T':
		case 't':
			myByte = NucleotideAlignmentConstants.T_ALLELE;
			break;
		case '+':
			myByte = NucleotideAlignmentConstants.INSERT_ALLELE;
		case '-':
			myByte = NucleotideAlignmentConstants.GAP_ALLELE;
		case 'N':
		case 'n':
			myByte = GenotypeTable.UNKNOWN_ALLELE_CHAR;
		default:
			myByte = NucleotideAlignmentConstants.UNDEFINED_ALLELE;
			break;
		}
		return myByte;
	}

}
