/*
 * MergeDuplicateSNPsPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.File;

import javax.swing.ImageIcon;
import org.apache.log4j.Logger;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.plugindef.AbstractPlugin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 * This class is intended to be run directly after TagsToSNPByAlignmentPlugin,
 * using the HapMap file from that step as input.
 *
 * It finds duplicate SNPs in the HapMap file, and merges them if they have the
 * same pair of alleles (not necessarily in the same maj/min order) and if there
 * mismatch rate is no greater than the threshold (-maxMisMat). If -callHets is
 * on, then genotypic disagreements will be called heterozygotes (otherwise set
 * to 'N' = default).
 *
 * By default, any remaining unmerged duplicate SNPs (but not indels) will be
 * deleted. They can be kept by invoking the -kpUnmergDups option.
 *
 * If the germplasm is not fully inbred, and still contains residual
 * heterozygosity (like the maize NAM or IBM populations do) then -callHets
 * should be on and -maxMisMat should be set fairly high (0.1 to 0.2, depending
 * on the amount of heterozygosity).
 *
 * @author jcg233
 */
public class MergeDuplicateSNPsPlugin extends AbstractPlugin {

    private static Logger myLogger = Logger.getLogger(MergeDuplicateSNPsPlugin.class);
    private static ArgsEngine myArgsEngine = null;
    private String suppliedInputFileName, suppliedOutputFileName, infile, outfile;
    private String snpLogFileName;
    private double maxMisMat = 0.05;
    private boolean usePedigree = false;
    HashMap<String, Double> taxaFs = null;
    boolean[] useTaxaForCompare = null;
    int nInbredTaxa = Integer.MIN_VALUE;
    private boolean callHets = false;  // true = when two genotypes disagree at a SNP, call it a heterozygote;  false = set to missing;
    private boolean kpUnmergDups = false;  // keep unmerged SNPs (not indels) in the data file
    int startChr = 1, endChr = 10;

    public MergeDuplicateSNPsPlugin() {
        super(null, false);
    }

    public MergeDuplicateSNPsPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
                "\n\nUsage is as follows:\n"
                + "-hmp           Input GBS genotype file (in HapMap format). Use a plus sign (+) as a wild card character to specify multiple chromosome numbers.\n"
                + "-o             Output HapMap file. Use a plus sign (+) as a wild card character to specify multiple chromosome numbers.\n"
                + "-misMat        Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: " + maxMisMat + ")\n"
                + "-p             Pedigree file containing full sample names (or expected names after merging) & expected inbreeding\n"
                + "                 coefficient (F) for each.  Only highly inbred taxa, with F >= 0.8 (e.g., S3 or more), will be used\n"
                + "                 to test if two duplicate SNPs agree with each other (default: use ALL taxa to compare duplicate SNPs)\n"
                + "-callHets      When two genotypes disagree at a SNP, call it a heterozygote (default: " + callHets + " = set to missing)\n"
                + "-kpUnmergDups  When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: " + kpUnmergDups + " = delete them)\n"
                + "-s             Start chromosome\n"
                + "-e             End chromosome\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-hmp", "--hmpFile", true);
            myArgsEngine.add("-o", "--outFile", true);
            myArgsEngine.add("-misMat", "--maxMismatchRate", true);
            myArgsEngine.add("-p", "--pedigree-file", true);
            myArgsEngine.add("-callHets", "--callHeterozygotes", false);
            myArgsEngine.add("-kpUnmergDups", "--keepUnmergedDuplicates", false);
            myArgsEngine.add("-s", "--startChromosome", true);
            myArgsEngine.add("-e", "--endChromosome", true);
            myArgsEngine.add("-snpLog", "", true);
        }
        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-hmp")) {
            suppliedInputFileName = myArgsEngine.getString("-hmp");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a HapMap file to filter.\n");
        }
        if (myArgsEngine.getBoolean("-o")) {
            suppliedOutputFileName = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file name.\n");
        }
        if (myArgsEngine.getBoolean("-misMat")) {
            maxMisMat = Double.parseDouble(myArgsEngine.getString("-misMat"));
        }
        if (myArgsEngine.getBoolean("-p")) {
            String pedigreeFileStr = myArgsEngine.getString("-p");
            File pedigreeFile = new File(pedigreeFileStr);
            if (!pedigreeFile.exists() || !pedigreeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the pedigree input file (-p option: " + pedigreeFileStr + ").");
            }
            taxaFs = TagsToSNPByAlignmentPlugin.readTaxaFsFromFile(pedigreeFile);
            if (taxaFs == null) {
                throw new IllegalArgumentException("Problem reading the pedigree file. Progam aborted.");
            }
            usePedigree = true;
        }
        if (myArgsEngine.getBoolean("-callHets")) {
            callHets = true;
        }
        if (myArgsEngine.getBoolean("-kpUnmergDups")) {
            kpUnmergDups = true;
        }
        if (myArgsEngine.getBoolean("-s")) {
            startChr = Integer.parseInt(myArgsEngine.getString("-s"));
        }
        if (myArgsEngine.getBoolean("-e")) {
            endChr = Integer.parseInt(myArgsEngine.getString("-e"));
        }
        if (endChr - startChr < 0) {
            printUsage();
            throw new IllegalArgumentException("Error: The start chromosome is higher than the end chromosome.");
        }
        if (myArgsEngine.getBoolean("-snpLog")) {
            snpLogFileName = myArgsEngine.getString("-snpLog");
        }
    }

    public DataSet performFunction(DataSet input) {
        for (int chr = startChr; chr <= endChr; chr++) {
            infile = suppliedInputFileName.replace("+", "" + chr);
            outfile = suppliedOutputFileName.replace("+", "" + chr);
            myLogger.info("Reading: " + infile);
            Alignment a;
            try {
                a = ImportUtils.readFromHapmap(infile, this);
                myLogger.info("Original Alignment  Taxa:" + a.getSequenceCount() + " Sites:" + a.getSiteCount());
                if (usePedigree && !maskNonInbredTaxa(a)) {
                    throw new IllegalArgumentException("Mismatch between taxa names in the input hapmap and pedigree files.");
                }
            } catch (Exception e) {
                myLogger.info("Could not read input hapmap file for chr" + chr + ":\n\t" + infile + "\n\te: " + e + "\n\tSkipping...");
                continue;
            }
            MutableNucleotideAlignment msa = MutableNucleotideAlignment.getInstance(a.getIdGroup(), a.getSiteCount());
            //MutableNucleotideAlignment msa = new MutableNucleotideAlignment(a.getIdGroup(), a.getSiteCount(), a.getLoci());
            ArrayList<Integer> samePosAL = new ArrayList<Integer>();
            Integer[] samePos = null;
            int currentPos = a.getPositionInLocus(0);
            for (int s = 0; s < a.getSiteCount(); s++) {  // must be sorted by position, as HapMap files typically are (ImportUtils.readFromHapmap() fails if they aren't)
                int newPos = a.getPositionInLocus(s);
                if (newPos == currentPos) {   // assumes that the strands are all '+' (in TagsToSNPByAlignmentPlugin(), - strand genos were complemented)
                    samePosAL.add(s);  // collect markers with the same position
                } else {
                    samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
                    if (samePosAL.size() > 1) {  // merge sets of 2 or more markers with the same position and alleles (alleles are not necessarily in the same order, maj/min)
                        processSNPsWithSamePosition(samePos, a, chr, currentPos, msa);
                    } else {  // site has a unique position: write its genos to the msa
                        byte[] genos = new byte[a.getSequenceCount()];
                        for (int t = 0; t < a.getSequenceCount(); ++t) {
                            genos[t] = a.getBase(t, samePos[0]);
                        }
                        addSiteToMutableAlignment(chr, currentPos, genos, msa);
                    }
                    // start a new collection of markers
                    samePosAL = new ArrayList<Integer>();
                    samePosAL.add(s);
                    currentPos = newPos;
                }
            }
            // finish last site or set of sites
            samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
            if (samePosAL.size() > 1) {
                processSNPsWithSamePosition(samePos, a, chr, currentPos, msa);
            } else {  // site has a unique position: write its genos to the msa
                byte[] genos = new byte[a.getSequenceCount()];
                for (int t = 0; t < a.getSequenceCount(); ++t) {
                    genos[t] = a.getBase(t, samePos[0]);
                }
                addSiteToMutableAlignment(chr, currentPos, genos, msa);
            }
            msa.clean();
            myLogger.info("Number of sites written after merging duplicate SNPs: " + msa.getSiteCount());
            if (!kpUnmergDups) {
                deleteRemainingDuplicates(msa);
            }
            ExportUtils.writeToHapmap(msa, false, outfile, '\t', this);
        }
        return null;
    }

    private void processSNPsWithSamePosition(Integer[] samePos, Alignment a, int chr, int currentPos, MutableNucleotideAlignment msa) {
        boolean[] finished = new boolean[samePos.length];
        for (int i = 0; i < finished.length; ++i) {
            finished[i] = false;   // indicates if the site has already been merged with a previous site OR written as is to msa
        }
        for (int s1 = 0; s1 < samePos.length - 1; ++s1) {
            if (finished[s1]) {
                continue;  // s1 has already been merged with a previous site
            }
            byte[] currMerge = new byte[a.getSequenceCount()];
            for (int t = 0; t < a.getSequenceCount(); ++t) {
                currMerge[t] = a.getBase(t, samePos[s1]); // set the current merger of genotypes to those for site s1
            }
            byte[] currAlleles = new byte[2];
            currAlleles[0] = a.getMajorAllele(samePos[s1].intValue());
            currAlleles[1] = a.getMinorAllele(samePos[s1].intValue());
            if (currAlleles[0] == NucleotideAlignmentConstants.GAP_ALLELE || currAlleles[1] == NucleotideAlignmentConstants.GAP_ALLELE) {
                addSiteToMutableAlignment(chr, currentPos, currMerge, msa);
                finished[s1] = true;
                continue;
            }
            Arrays.sort(currAlleles);
            for (int s2 = s1 + 1; s2 < samePos.length; ++s2) {
                if (finished[s2]) {
                    continue;  // s2 has already been merged with a previous site (perhaps with different alleles)
                }
                byte[] newAlleles = new byte[2];
                newAlleles[0] = a.getMajorAllele(samePos[s2].intValue());
                newAlleles[1] = a.getMinorAllele(samePos[s2].intValue());
                Arrays.sort(newAlleles);
                if (newAlleles[0] == NucleotideAlignmentConstants.GAP_ALLELE || newAlleles[1] == NucleotideAlignmentConstants.GAP_ALLELE) {
                    continue;
                }
                // Check if the alleles match.  If they do, merge the genos, provided that the number of genotypic mismatches is below threshold
                if (Arrays.equals(currAlleles, newAlleles)) {
                    int nMismatch = 0;
                    int nCompare = 0;
                    byte[] possibleMerge = new byte[a.getSequenceCount()];
                    for (int t = 0; t < a.getSequenceCount(); ++t) {
                        byte geno2 = a.getBase(t, samePos[s2]);
                        if (currMerge[t] != Alignment.UNKNOWN_DIPLOID_ALLELE && geno2 != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            if (!usePedigree || useTaxaForCompare[t]) {
                                ++nCompare;
                            }
                            if (!AlignmentUtils.isEqual(currMerge[t], geno2)) {
                                if (!usePedigree || useTaxaForCompare[t]) {
                                    ++nMismatch;
                                }
                                try {
                                    possibleMerge[t] = callHets ? resolveHet(currMerge[t], geno2) : Alignment.UNKNOWN_DIPLOID_ALLELE;
                                } catch (Exception e) {
                                    myLogger.warn(
                                            "Invalid genotypes (" + a.getBaseAsString(t, samePos[s1]) + " and " + a.getBaseAsString(t, samePos[s2]) + ") at position:" + currentPos + " taxon:" + a.getTaxaName(t));
                                }
                            } else {
                                possibleMerge[t] = currMerge[t];
                            }
                        } else if (currMerge[t] == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            possibleMerge[t] = geno2;
                        } else if (geno2 == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            possibleMerge[t] = currMerge[t];
                        }
                    }
                    if ((nCompare == 0) || ((double) nMismatch / nCompare <= maxMisMat)) {
                        for (int t = 0; t < a.getSequenceCount(); ++t) {
                            currMerge[t] = possibleMerge[t];
                        }
                        finished[s2] = true;
                    }
                }
            }
            // Finished comparing s1 to all other sites
            // Write the current "merged" genos to the msa (if all comparisons disagreed, these will be the original s1 genos)
            addSiteToMutableAlignment(chr, currentPos, currMerge, msa);
            finished[s1] = true;
        }
        // make sure that everything in the current collection (samePos) got finished
        //  (for example, the final site may not have been merged with anything, or might have had gap as an allele)
        for (int site = 0; site < finished.length; ++site) {
            if (!finished[site]) {
                byte[] genos = new byte[a.getSequenceCount()];
                for (int t = 0; t < a.getSequenceCount(); ++t) {
                    genos[t] = a.getBase(t, samePos[site]);
                }
                addSiteToMutableAlignment(chr, currentPos, genos, msa);
            }
        }
    }

    private void addSiteToMutableAlignment(int chromosome, int position, byte[] genos, MutableNucleotideAlignment theMSA) {
        int currSite = theMSA.getSiteCount();
        //int currSite = theMSA.getNextFreeSite();
        theMSA.addSite(currSite);
        theMSA.setLocusOfSite(currSite, new Locus(String.valueOf(chromosome), String.valueOf(chromosome), -1, -1, null, null));
        theMSA.setPositionOfSite(currSite, position);
        //theMSA.setStrandOfSite(currSite, (byte) '+');
        for (int tx = 0; tx < genos.length; tx++) {
            theMSA.setBase(tx, currSite, genos[tx]);
        }
    }

    private void deleteRemainingDuplicates(MutableNucleotideAlignment theMSA) {
        SNPLogging snpLogging = new SNPLogging(snpLogFileName);
        ArrayList<Integer> samePosAL = new ArrayList<Integer>();
        int currentPos = theMSA.getPositionInLocus(0);
        for (int s = 0; s < theMSA.getSiteCount(); s++) {
            int newPos = theMSA.getPositionInLocus(s);
            if (newPos == currentPos) {
                samePosAL.add(s);  // collect markers with the same position
            } else {
                if (samePosAL.size() > 1) {
                    Integer[] samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
                    for (int i = 0; i < samePos.length; ++i) {
                        if (theMSA.getMajorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE
                                && theMSA.getMinorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
                            snpLogging.writeEntry(theMSA, samePos[i], null, null, this.getClass(), null, null, null, null);
                            theMSA.clearSiteForRemoval(samePos[i]);
                        }
                    }
                }
                // start a new collection of markers
                samePosAL = new ArrayList<Integer>();
                samePosAL.add(s);
                currentPos = newPos;
            }
        }
        // finish last site or set of sites
        if (samePosAL.size() > 1) {
            Integer[] samePos = samePosAL.toArray(new Integer[samePosAL.size()]);
            for (int i = 0; i < samePos.length; ++i) {
                if (theMSA.getMajorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE
                        && theMSA.getMinorAllele(samePos[i].intValue()) != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
                    snpLogging.writeEntry(theMSA, samePos[i], null, null, this.getClass(), null, null, null, null);
                    theMSA.clearSiteForRemoval(samePos[i]);
                }
            }
        }
        snpLogging.close();
        theMSA.clean();
        myLogger.info("Number of sites written after deleting any remaining, unmerged duplicate SNPs: " + theMSA.getSiteCount());
    }

    private boolean maskNonInbredTaxa(Alignment a) {
        useTaxaForCompare = new boolean[a.getSequenceCount()];  // initialized to false
        nInbredTaxa = 0;
        try {
            for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
                if (taxaFs.containsKey(a.getFullTaxaName(taxon))) {
                    if (taxaFs.get(a.getFullTaxaName(taxon)) >= 0.8) {
                        useTaxaForCompare[taxon] = true;
                        nInbredTaxa++;
                    }
                } else {
                    throw new Exception("Taxon " + a.getFullTaxaName(taxon) + " not found in the pedigree file");
                }
            }
            myLogger.info(nInbredTaxa + " highly inbred taxa (with an expected F >= 0.8) were found in the input hapmap file (according to the pedigree file)");
            return true;
        } catch (Exception e) {
            myLogger.error("Mismatch between taxa names in the input hapmap file and the pedigree file e=" + e);
            e.printStackTrace();
            return false;
        }
    }

    private static byte resolveHet(byte geno1, byte geno2) {

        byte[] result = new byte[2];
        result[0] = (byte) (geno1 >>> 4);
        byte temp = (byte) (geno1 & 0xf);
        int count = 1;
        if (temp != result[0]) {
            result[count++] = temp;
        }

        temp = (byte) (geno2 >>> 4);
        if (temp == result[0]) {
            // do nothing
        } else if (count == 1) {
            result[count++] = temp;
        } else if (temp != result[1]) {
            throw new IllegalStateException();
        }

        temp = (byte) (geno2 & 0xf);
        if (temp == result[0]) {
            // do nothing
        } else if (count == 1) {
            result[count++] = temp;
        } else if (temp != result[1]) {
            throw new IllegalStateException();
        }

        return (byte) ((result[0] << 4) | (result[1]));

    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
