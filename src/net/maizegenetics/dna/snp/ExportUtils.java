/*
 * ExportUtils
 */
package net.maizegenetics.dna.snp;

import com.google.common.base.Joiner;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Ints;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

/**
 * Exports Genotype Tables to various file formats.
 *
 * @author Jon Zhang
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class ExportUtils {

    private static final Logger myLogger = Logger.getLogger(ExportUtils.class);
    private static FormattedOutput format = FormattedOutput.getInstance();

    private ExportUtils() {
        // Utility Class - do not instantiate.
    }

    public static String writeGenotypeHDF5(GenotypeTable a, String newHDF5file) {
        return writeGenotypeHDF5(a, newHDF5file, null, true);
    }

    public static String writeGenotypeHDF5(GenotypeTable a, String newHDF5file, boolean keepDepth) {
        return writeGenotypeHDF5(a, newHDF5file, null, keepDepth);
    }

    /**
     * Exports a alignment into the Byte HDF5 format.
     *
     * @param a alignment to be exported
     * @param newHDF5file filename for the new file (should end with "hmp.h5")
     * @param exportTaxa subset of taxa (if null exports ALL taxa)
     * @return
     */
    public static String writeGenotypeHDF5(GenotypeTable a, String newHDF5file, TaxaList exportTaxa, boolean keepDepth) {
        GenotypeTableBuilder aB = GenotypeTableBuilder.getTaxaIncremental(a.positions(), newHDF5file);
        if ((exportTaxa != null) && (exportTaxa.numberOfTaxa() == 0)) {
            aB.build();
            return newHDF5file;
        }
        for (int t = 0; t < a.numberOfTaxa(); t++) {
            if ((exportTaxa != null) && (!exportTaxa.contains(a.taxa().get(t)))) {
                continue;  //taxon not in export list
            }
            byte[] bases = a.genotypeAllSites(t);
            if (a.hasDepth() == false || keepDepth == false) {
                aB.addTaxon(a.taxa().get(t), bases, null);
            } else {
                aB.addTaxon(a.taxa().get(t), bases, a.depth().depthAllSitesByte(t));
            }
        }
        aB.build();
        return newHDF5file;
    }

    /**
     * Write a GenotypeTable to HapMap format with standard settings - unphased
     * single character, tab delimiter, and no progress tracking.
     *
     * @param alignment genotype table
     * @param filename outfile name (will add ".hmp.txt" if needed)
     * @return name of the outfile with the appropriate suffix
     */
    public static String writeToHapmap(GenotypeTable alignment, String filename) {
        return writeToHapmap(alignment, false, filename, '\t', null);
    }

    /**
     * Write a GenotypeTable to HapMap format.
     *
     * @param alignment genotype table
     * @param diploid true uses phased two letter encoding, false one letter
     * unphased
     * @param filename outfile name (will add ".hmp.txt" if needed)
     * @param delimChar delimiter character normally tab
     * @param listener progress listener, (null if unneeded)
     * @return name of the outfile with the appropriate suffix
     */
    public static String writeToHapmap(GenotypeTable alignment, boolean diploid, String filename, char delimChar, ProgressListener listener) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(filename, ".hmp.txt", new String[]{".hmp.txt", ".hmp.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
            if (true) {
                for (Taxon taxon : alignment.taxa()) {
                    if (taxon.getAllAnnotationEntries().length == 0) {
                        continue;
                    }
                    bw.write("##SAMPLE=" + taxon.toStringWithVCFAnnotation() + "\n");
                }
            }
            bw.write(Joiner.on(delimChar).join("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID",
                    "assayLSID", "panelLSID", "QCcode"));
            bw.write(delimChar);
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String sequenceID = alignment.taxaName(taxa).trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                bw.write(alignment.siteName(site));
                bw.write(delimChar);
                int[][] sortedAlleles = alignment.allelesSortedByFrequency(site); // which alleles are actually present among the genotypes
                int numAlleles = sortedAlleles[0].length;
                if (numAlleles == 0) {
                    bw.write("NA"); //if data does not exist
                } else if (numAlleles == 1) {
                    bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                } else {
                    bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                    for (int allele = 1; allele < sortedAlleles[0].length; allele++) {
                        if (sortedAlleles[0][allele] != GenotypeTable.UNKNOWN_ALLELE) {
                            bw.write('/');
                            bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][allele]));  // will write out a third allele if it exists
                        }
                    }
                }
                bw.write(delimChar);
                bw.write(Joiner.on(delimChar).join(alignment.chromosomeName(site), String.valueOf(alignment.chromosomalPosition(site)),
                        "+", "NA", "NA", "NA", "NA", "NA", "NA"));
                bw.write(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
                        String baseIUPAC = null;
                        try {
                            baseIUPAC = alignment.genotypeAsString(taxa, site);
                        } catch (Exception e) {
                            String[] b = alignment.genotypeAsStringArray(taxa, site);
                            bw.close();
                            myLogger.debug(e.getMessage(), e);
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                        }
                        if ((baseIUPAC == null) || baseIUPAC.equals("?")) {
                            String[] b = alignment.genotypeAsStringArray(taxa, site);
                            bw.close();
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                        }
                        bw.write(baseIUPAC);
                    } else {
                        String[] b = alignment.genotypeAsStringArray(taxa, site);
                        if (b.length == 1) {
                            bw.write(b[0]);
                            bw.write(b[0]);
                        } else {
                            bw.write(b[0]);
                            bw.write(b[1]);
                        }
                    }
                    if (taxa != (numTaxa - 1)) {
                        bw.write(delimChar);
                    }
                }
                bw.write("\n");

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }
            return fullFileName;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("Error writing Hapmap file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Writes given alignment to a VCF file
     *
     * @param gt
     * @param filename
     * @return
     */
    public static String writeToVCF(GenotypeTable gt, String filename, boolean keepDepth) {
        final char delimChar = '\t';
        boolean hasDepth = gt.hasDepth() && keepDepth;
        try {

            filename = Utils.addSuffixIfNeeded(filename, ".vcf", new String[]{".vcf", ".vcf.gz"});
            BufferedWriter bw = Utils.getBufferedWriter(filename);
            bw.write("##fileformat=VCFv4.0");
            bw.newLine();
            if (!gt.hasReference()) {
                bw.write("##Tassel=<ID=GenotypeTable,Version=5,Description=\"Reference allele is not known. The major allele was used as reference allele\">");
                bw.newLine();
            }
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">");
            bw.newLine();
            bw.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
            bw.newLine();
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
            bw.newLine();
            bw.write("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">");
            bw.newLine();
            writeVCFSampleAnnotationToWriter(gt, bw);
            bw.write("#CHROM" + delimChar + "POS" + delimChar + "ID" + delimChar + "REF" + delimChar + "ALT" + delimChar + "QUAL" + delimChar + "FILTER" + delimChar + "INFO" + delimChar + "FORMAT");
            for (int taxa = 0; taxa < gt.numberOfTaxa(); taxa++) {
                String taxonName = gt.taxaName(taxa).trim();
                bw.write(delimChar + taxonName);
            }
            bw.newLine();

            int noAlleles = 0;
            for (int site = 0; site < gt.numberOfSites(); site++) {
                Position p = gt.positions().get(site);
                byte refAllele = p.getAllele(WHICH_ALLELE.Reference);
                int[] sortedAlleles = gt.allelesSortedByFrequency(site)[0]; // which alleles are actually present among the genotypes
                int indexOfRefAllele = Ints.indexOf(sortedAlleles, refAllele);
                if (indexOfRefAllele < 0) {
                    indexOfRefAllele = 0;
                }
                if (indexOfRefAllele != 0) {
                    int t = sortedAlleles[0];
                    sortedAlleles[0] = sortedAlleles[indexOfRefAllele];
                    sortedAlleles[indexOfRefAllele] = t;
                }

                int nAlleles = sortedAlleles.length;
                String[] alleleRedirect = new String[16];
                Arrays.fill(alleleRedirect, ".");
                for (int i = 0; i < sortedAlleles.length; i++) {
                    alleleRedirect[sortedAlleles[i]] = "" + i;
                }

                bw.write(gt.chromosomeName(site)); // chromosome
                bw.write(delimChar);
                bw.write(gt.chromosomalPosition(site) + ""); // position
                bw.write(delimChar);
                bw.write(gt.siteName(site)); // site name
                bw.write(delimChar);
                if (nAlleles == 0) {                                                  //used to be ==0
                    //System.out.println("A0:"+gt.chromosomeName(site)+":"+gt.chromosomalPosition(site));
                    noAlleles++;
                    bw.write(".\t.\t.\tPASS\t.\tGT");
                    for (int taxa = 0; taxa < gt.numberOfTaxa(); taxa++) {
                        bw.write("\t./.");
                    }
                    bw.newLine();
                    continue;
                }
                bw.write(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[0])); // ref allele
                bw.write(delimChar);

                StringBuilder altAllelesBuilder = new StringBuilder("");
                for (int aa = 1; aa < sortedAlleles.length; aa++) {
                    altAllelesBuilder.append(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[aa]) + ",");
                }
                if (altAllelesBuilder.length() == 0) {
                    altAllelesBuilder.append(".");
                } else {
                    altAllelesBuilder.deleteCharAt(altAllelesBuilder.length() - 1);
                }
                bw.write(altAllelesBuilder.toString()); // alt alleles
                bw.write(delimChar);

                bw.write("."); // qual score
                bw.write(delimChar);

                bw.write("PASS"); // filter
                bw.write(delimChar);

                if (hasDepth) {
                    bw.write("DP=" + gt.depth().depthForSite(site)); // DP
                } else {
                    bw.write("."); // DP
                }
                bw.write(delimChar);

                if (hasDepth) {
                    bw.write("GT:AD:DP:GQ:PL");
                } else {
                    bw.write("GT");
                }
                for (int taxa = 0; taxa < gt.numberOfTaxa(); taxa++) {
                    bw.write(delimChar);
                    // GT = genotype
                    byte[] values = gt.genotypeArray(taxa, site);
                    bw.write(alleleRedirect[values[0]]);
                    bw.write("/");
                    bw.write(alleleRedirect[values[1]]);
                    if (!(hasDepth)) {
                        continue;
                    }
                    bw.write(":");

                    // AD
                    int[] siteAlleleDepths = gt.depthForAlleles(taxa, site);

                    int siteTotalDepth = 0;
                    for (int ss = 0; ss < sortedAlleles.length; ss++) {
                        bw.write("" + siteAlleleDepths[sortedAlleles[ss]]);
                        if (ss < sortedAlleles.length - 1) {
                            bw.write(',');
                        }
                        siteTotalDepth += siteAlleleDepths[sortedAlleles[ss]];
                    }
                    bw.write(":");
                    // DP
                    bw.write(siteTotalDepth + "");
                    bw.write(":");
                    int[] scores = new int[]{-1, -1, -1, -1};
                    if (values[0] != GenotypeTable.UNKNOWN_ALLELE) {
                        int altDepth = (sortedAlleles.length < 2) ? 0 : siteAlleleDepths[sortedAlleles[1]];
                        scores = VCFUtil.getScore(siteAlleleDepths[sortedAlleles[0]], altDepth);
                    }  //this is not correct
                    //  else {scores= VCFUtil.getScore(siteAlleleDepths[values[0]], siteAlleleDepths[values[1]]);} //this is not correct
                    // GQ
                    bw.write(scores[3] + "");
                    bw.write(":");
                    // PL
                    bw.write(scores[0] + "," + scores[1] + "," + scores[2]);
//                    } else {
//                        bw.write(".:.,.");
//                    }
                }
                bw.newLine();
            }
            if (noAlleles > 0) {
                myLogger.warn("Warning: " + noAlleles + " sites have no alleles.");
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing VCF file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        }
        return filename;
    }

    private static void writeVCFSampleAnnotationToWriter(GenotypeTable gt, BufferedWriter bw) throws IOException {
        for (Taxon taxon : gt.taxa()) {
            Multimap annoMap = taxon.getAnnotationAsMap();
            if (annoMap.size() == 0) {
                continue;
            }
            String annoString = Joiner.on(',').withKeyValueSeparator("=").join(annoMap.entries());
            bw.write("##SAMPLE=<ID=" + taxon.getName() + "," + annoString + ">");
            bw.newLine();
        }
    }

    /**
     * Writes given set of alignments to a set of Plink files
     *
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static String writeToPlink(GenotypeTable alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter MAPbw = null;
        BufferedWriter PEDbw = null;
        String mapFileName = Utils.addSuffixIfNeeded(filename, ".plk.map");
        String pedFileName = Utils.addSuffixIfNeeded(filename, ".plk.ped");
        try {
            MAPbw = new BufferedWriter(new FileWriter(mapFileName), 1000000);
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.chromosomeName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(alignment.siteName(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write("-9"); // genetic distance unavailable
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.chromosomalPosition(site))); // position
                MAPbw.write("\n");
            }
            MAPbw.close();

            PEDbw = new BufferedWriter(new FileWriter(pedFileName), 1000000);
            // Compiled : Pattern
            Pattern splitter = Pattern.compile(":");
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] name = splitter.split(alignment.taxaName(taxa).trim());
                if (name.length != 1) {
                    PEDbw.write(name[1]); // namelvl 1 if is available
                } else {
                    PEDbw.write("-9");
                }
                PEDbw.write(delimChar);
                PEDbw.write(alignment.taxaName(taxa).trim()); // namelvl 0
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // paternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // maternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // gender is both
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // phenotype unavailable, might have to change to "-9" for missing affection status
                PEDbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    String[] b = getSNPValueForPlink(alignment.genotypeAsStringArray(taxa, site));
                    PEDbw.write(b[0]);
                    PEDbw.write(delimChar);
                    PEDbw.write(b[b.length - 1]);
                    if (site != numSites - 1) {
                        PEDbw.write(delimChar);
                    }
                }
                PEDbw.write("\n");
            }
            PEDbw.close();
            return mapFileName + " and " + pedFileName;
        } catch (Exception e) {
            myLogger.error("Error writing Plink files: " + mapFileName + " and " + pedFileName + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Plink files: " + mapFileName + " and " + pedFileName + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                PEDbw.close();
            } catch (Exception e) {
                // do nothing
            }
            try {
                MAPbw.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    private static String[] getSNPValueForPlink(String[] base) {
        for (int i = 0; i < base.length; i++) {
            if (base[i].equals("N")) {
                base[i] = "0";
            } else if (base[i].equals("0")) {
                base[i] = "D";
            }
        }
        return base;
    }

    public static String saveDelimitedAlignment(GenotypeTable theAlignment, String delimit, String saveFile) {

        if ((saveFile == null) || (saveFile.length() == 0)) {
            return null;
        }
        saveFile = Utils.addSuffixIfNeeded(saveFile, ".txt");
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(new File(saveFile));
            bw = new BufferedWriter(fw);

            bw.write("Taxa");
            int numSites = theAlignment.numberOfSites();
            for (int j = 0; j < numSites; j++) {
                bw.write(delimit);
                bw.write(String.valueOf(theAlignment.chromosomalPosition(j)));
            }
            bw.write("\n");

            for (int r = 0, n = theAlignment.numberOfTaxa(); r < n; r++) {
                bw.write(theAlignment.taxaName(r));
                for (int i = 0; i < numSites; i++) {
                    bw.write(delimit);
                    bw.write(theAlignment.genotypeAsString(r, i));
                }
                bw.write("\n");
            }

            return saveFile;

        } catch (Exception e) {
            myLogger.error("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    /**
     * print alignment (in PHYLIP SEQUENTIAL format)
     */
    public static void printSequential(GenotypeTable a, PrintWriter out) {
        // PHYLIP header line
        out.println("  " + a.numberOfTaxa() + " " + a.numberOfSites() + "  S");

        // Print sequences
        for (int s = 0; s < a.numberOfTaxa(); s++) {
            int n = 0;
            while (n < a.numberOfSites()) {
                if (n == 0) {
                    format.displayLabel(out, a.taxaName(s), 10);
                    out.print("     ");
                } else {
                    out.print("               ");
                }
                printNextSites(a, out, false, s, n, 50);
                out.println();
                n += 50;
            }
        }
    }

    /**
     * print alignment (in PHYLIP 3.4 INTERLEAVED format)
     */
    public static void printInterleaved(GenotypeTable a, PrintWriter out) {
        int n = 0;

        // PHYLIP header line
        out.println("  " + a.numberOfTaxa() + " " + a.numberOfSites());

        // Print sequences
        while (n < a.numberOfSites()) {
            for (int s = 0; s < a.numberOfTaxa(); s++) {
                if (n == 0) {
                    format.displayLabel(out, a.taxaName(s), 10);
                    out.print("     ");
                } else {
                    out.print("               ");
                }
                printNextSites(a, out, true, s, n, 50);
                out.println();
            }
            out.println();
            n += 50;
        }
    }

    /**
     * Print alignment (in CLUSTAL W format)
     */
    public static void printCLUSTALW(GenotypeTable a, PrintWriter out) {
        int n = 0;

        // CLUSTAL W header line
        out.println("CLUSTAL W multiple sequence alignment");
        out.println();

        // Print sequences
        while (n < a.numberOfSites()) {
            out.println();
            for (int s = 0; s < a.numberOfTaxa(); s++) {
                format.displayLabel(out, a.taxaName(s), 10);
                out.print("     ");

                printNextSites(a, out, false, s, n, 50);
                out.println();
            }
            // Blanks in status line are necessary for some parsers)
            out.println("               ");
            n += 50;
        }
    }

    private static void printNextSites(GenotypeTable a, PrintWriter out, boolean chunked, int seq, int start, int num) {
        // Print next num characters
        for (int i = 0; (i < num) && (start + i < a.numberOfSites()); i++) {
            // Chunks of 10 characters
            if (i % 10 == 0 && i != 0 && chunked) {
                out.print(' ');
            }
            out.print(a.genotypeAsString(seq, start + i));
        }
    }

    public static String writeAlignmentToSerialGZ(GenotypeTable sba, String outFile) {

        long time = System.currentTimeMillis();

        File theFile = null;
        FileOutputStream fos = null;
        GZIPOutputStream gz = null;
        ObjectOutputStream oos = null;
        try {
            theFile = new File(Utils.addSuffixIfNeeded(outFile, ".serial.gz"));
            fos = new FileOutputStream(theFile);
            gz = new GZIPOutputStream(fos);
            oos = new ObjectOutputStream(gz);
            oos.writeObject(sba);
            return theFile.getName();
        } catch (Exception e) {
            e.printStackTrace();
            myLogger.error("Error writing Serial GZ: " + theFile.getName() + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Serial GZ: " + theFile.getName() + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                oos.flush();
                oos.close();
                gz.close();
                fos.close();
            } catch (Exception e) {
                // do nothing
            }
            myLogger.info("writeAlignmentToSerialGZ: " + theFile.toString() + "  Time: " + (System.currentTimeMillis() - time));
        }

    }
}
