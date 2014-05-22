/*
 * AddReferenceAlleleToHDF5Plugin
 */
package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;
import net.maizegenetics.util.ArgsEngine;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Utils;
import java.util.List;

/**
 *
 * @author Jeff Glaubitz (jcg233@cornell.edu)
 */
public class AddReferenceAlleleToHDF5Plugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(AddReferenceAlleleToHDF5Plugin.class);
    private ArgsEngine myArgsEngine = null;
    private String inHDF5FileName = null;
    private String refGenomeFileStr = null;
    private String genomeVersion = null;
    private BufferedReader refReader = null;
    private PositionListBuilder newPosListBuilder = null;
    private int currChr = Integer.MIN_VALUE;
    private int currPos = Integer.MIN_VALUE;
    private String contextSeq;
    private final boolean writePositions = false;

    public AddReferenceAlleleToHDF5Plugin() {
        super(null, false);
    }

    public AddReferenceAlleleToHDF5Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
            "\n\n"
            +"The options for the TASSEL5 AddReferenceAlleleToHDF5Plugin are as follows:\n"
            +"  -i    Input HDF5 genotype (*.h5) file to be annotated with the reference allele.\n"
            +"  -ref  Path to reference genome in fasta format.\n"
            +"  -ver  Genome version.\n"
            +"\n\n"
        );
    }

    @Override
    public void setParameters(String[] args) {
        if (args == null || args.length == 0) {
            printUsage();
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-file", true);
            myArgsEngine.add("-ref", "--referenceGenome", true);
            myArgsEngine.add("-ver", "--genome-version", true);
        }
        myArgsEngine.parse(args);
        inHDF5FileName = myArgsEngine.getString("-i");
        File inFile = null;
        if (inHDF5FileName != null) {
            inFile = new File(inHDF5FileName);
            if (!inFile.exists()) {
                printUsage();
                try {Thread.sleep(500);} catch(InterruptedException e) {}
                String errMessage = "\nThe input HDF5 genotype (*.h5) file name you supplied does not exist:\n"+inHDF5FileName+"\n\n";
                myLogger.error(errMessage);
                try {Thread.sleep(500);} catch(InterruptedException e) {}
                throw new IllegalArgumentException(errMessage);
            }
        } else {
            printUsage();
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            String errMessage = "\nPlease specify an input HDF5 genotype (*.h5) file to be annotated with the reference allele (-i option)\n\n";
            myLogger.error(errMessage);
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            throw new IllegalArgumentException(errMessage);
        }
        if (myArgsEngine.getBoolean("-ref")) {
            refGenomeFileStr = myArgsEngine.getString("-ref");
            File refGenomeFile = new File(refGenomeFileStr);
            if (!refGenomeFile.exists() || !refGenomeFile.isFile()) {
                printUsage();
                try {Thread.sleep(500);} catch(InterruptedException e) {}
                String errMessage = "\nCan't find the reference genome fasta file specified by the -ref option:\n  "+refGenomeFileStr+"\n\n";
                myLogger.error(errMessage);
                try {Thread.sleep(500);} catch(InterruptedException e) {}
                throw new IllegalArgumentException();
            }
            refReader = Utils.getBufferedReader(refGenomeFileStr);
        } else {
            printUsage();
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            String errMessage = "\nPlease specify a reference genome fasta file (-ref option)\n\n";
            myLogger.error(errMessage);
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            throw new IllegalArgumentException(errMessage);
        }
        if (myArgsEngine.getBoolean("-ver")) {
            genomeVersion = myArgsEngine.getString("-ver");
        } else {
            printUsage();
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            String errMessage = "\nPlease specify the name of the reference genome version (e.g. -ver \"B73 RefGen_v2\")\n\n";
            myLogger.error(errMessage);
            try {Thread.sleep(500);} catch(InterruptedException e) {}
            throw new IllegalArgumentException(errMessage);
        }
    }

    @Override
    public DataSet performFunction(DataSet input) {
        String message = addRefAlleleToHDF5GenoTable();
        if(message != null) {
            myLogger.error(message);
            try {Thread.sleep(500);} catch(Exception e) {}
            throw new IllegalStateException(message);
        }
        try {refReader.close();} catch (IOException e) {}
        fireProgress(100);
        return null;
    }

    private String addRefAlleleToHDF5GenoTable() {
        String message = populatePositionsWithRefAllele();
        if (message != null) return message;
        PositionList newPos = newPosListBuilder.build();
        if (writePositions) {
            String genomeVer = newPos.hasReference() ? newPos.genomeVersion() : "unkown";
            System.out.println("\nNew genome version: "+genomeVer+"\n");
        }
        String outHDF5FileName = inHDF5FileName.replaceFirst("\\.h5$", "_withRef.h5");
        GenotypeTableBuilder newGenos = GenotypeTableBuilder.getTaxaIncremental(newPos, outHDF5FileName);
        IHDF5Reader h5Reader = HDF5Factory.open(inHDF5FileName);
        List<String> taxaNames = HDF5Utils.getAllTaxaNames(h5Reader);
        for (String taxonName : taxaNames) {
            Taxon taxon = HDF5Utils.getTaxon(h5Reader, taxonName);
            byte[] genos = HDF5Utils.getHDF5GenotypesCalls(h5Reader, taxonName);
            byte[][] depth = HDF5Utils.getHDF5GenotypesDepth(h5Reader, taxonName);
            newGenos.addTaxon(taxon, genos, depth);
        }
        newGenos.build();
        System.out.println("\n\nFinished adding reference alleles to file:");
        System.out.println("  "+outHDF5FileName+"\n\n");
        return null;
    }
    
    private String populatePositionsWithRefAllele() {
        IHDF5Writer h5Writer = HDF5Factory.open(inHDF5FileName);
        PositionList oldPosList = PositionListBuilder.getInstance(h5Writer);
//        h5Writer.close();
        newPosListBuilder = new PositionListBuilder();
        if (writePositions) {
            String genomeVer = oldPosList.hasReference() ? oldPosList.genomeVersion() : "unkown";
            System.out.println("\ncurrent genome version: "+genomeVer+"\n");
            System.out.println("SNPID\tchr\tpos\tstr\tmaj\tmin\tref\tmaf\tcov\tcontext");
        }
        for (Position oldPos : oldPosList) {
            int chr = oldPos.getChromosome().getChromosomeNumber();
            int pos = oldPos.getPosition();
            byte strand = oldPos.getStrand();
            byte refAllele = retrieveRefAllele(chr, pos, strand);
            if (refAllele == Byte.MIN_VALUE) {
                return "\nCould not find position "+pos+" on chromosome "+chr+" in the reference genome fasta file.\n\n\n";
            }
            Position newPos = new GeneralPosition.Builder(oldPos) 
//                // previous version only copied chr,pos,strand,CM,SNPID,isNucleotide,isIndel from oldPos
//                .maf(oldPos.getGlobalMAF())
//                .siteCoverage(oldPos.getGlobalSiteCoverage())
//                .allele(Position.WHICH_ALLELE.GlobalMajor, oldPos.getAllele(Position.WHICH_ALLELE.GlobalMajor))
//                .allele(Position.WHICH_ALLELE.GlobalMinor, oldPos.getAllele(Position.WHICH_ALLELE.GlobalMinor))
//                .allele(Position.Allele.ANC, oldPos.getAllele(Position.Allele.ANC))
//                .allele(Position.Allele.HIDEP, oldPos.getAllele(Position.Allele.HIDEP))
                .allele(WHICH_ALLELE.Reference, refAllele)
                .build();
            if (writePositions) writePosition(newPos, contextSeq);
            newPosListBuilder.add(newPos);
        }
        newPosListBuilder.genomeVersion(genomeVersion);
        return null;
    }
    
    private byte retrieveRefAllele(int chr, int pos, int strand) {
        findChrInRefGenomeFile(chr);
        char currChar = findPositionInRefGenomeFile(pos);
        if (currPos == pos) {
            byte refAllele = NucleotideAlignmentConstants.getNucleotideAlleleByte(currChar);
            if (strand == -1) refAllele = NucleotideAlignmentConstants.getNucleotideComplement(refAllele);
            return refAllele;
        } else {
            System.out.println("currPos:"+currPos);
            return Byte.MIN_VALUE;
        }
    }
    
    private void findChrInRefGenomeFile(int chr) {
        String temp = "Nothing has been read from the reference genome fasta file yet";
        try {
            while (refReader.ready() && currChr < chr) {
                temp = refReader.readLine().trim();
                if (temp.startsWith(">")) {
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    try {
                        currChr = Integer.parseInt(chrS);
                    } catch (NumberFormatException e) {
                        myLogger.error("\n\nAddReferenceAlleleToHDF5Plugin detected a non-numeric chromosome name in the reference genome sequence fasta file: " + chrS
                                + "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
                                + "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n");
                        System.exit(1);
                    }
                    if (!writePositions) myLogger.info("\nCurrently reading chromosome "+currChr+" from reference genome fasa file\n\n");
                }
                currPos = 0;
            }
        } catch (IOException e) {
            myLogger.error("Exception caught while reading the reference genome fasta file:\n  "+e+"\nLast line read:\n  "+temp);
            try {Thread.sleep(500);} catch(InterruptedException iE) {}
            e.printStackTrace();
            System.exit(1);
        }
        if (currChr != chr) {
            myLogger.error("\nCould not find chromosome "+chr+" in the reference genome fasta file.\nMake sure that the chromosomes are in numerical order in that file\n\n\n");
            try {Thread.sleep(500);} catch(InterruptedException iE) {}
            System.exit(1);
        }
    }
    
    private char findPositionInRefGenomeFile(int pos) {
        char currChar = Character.MAX_VALUE;
        contextSeq = "";
        try {
            while (currPos < pos) {
                int intChar = refReader.read();
                if (intChar == -1) {
                    currPos = Integer.MAX_VALUE;
                    return Character.MAX_VALUE;
                }
                currChar = (char) intChar;
                if (!Character.isWhitespace(currChar)) {
                    currPos++;
                    if (pos - currPos < 60) {
                        contextSeq += currChar;
                    }
                }
            }
        } catch (IOException e) {
            myLogger.error("\n\nError reading reference genome file:\n  "+e+"\n\n");
            System.exit(1);
        }
        return currChar;
    }
    
    private void writePosition(Position pos, String contextSeq) {
        System.out.println(
            pos.getSNPID()+
            "\t"+pos.getChromosome().getChromosomeNumber()+
            "\t"+pos.getPosition()+
            "\t"+pos.getStrand()+
            "\t"+pos.getAllele(WHICH_ALLELE.GlobalMajor)+
            "\t"+pos.getAllele(WHICH_ALLELE.GlobalMinor)+
            "\t"+pos.getAllele(WHICH_ALLELE.Reference)+
            "\t"+pos.getGlobalMAF()+
            "\t"+pos.getGlobalSiteCoverage()+
            "\t"+contextSeq
        );
    }

    @Override
    public ImageIcon getIcon() {
//        URL imageURL = SeparatePlugin.class.getResource("/net/maizegenetics/analysis/images/Merge.gif");
//        if (imageURL == null) {
            return null;
//        } else {
//            return new ImageIcon(imageURL);
//        }
    }

    @Override
    public String getButtonName() {
        return "Add reference allele";
    }

    @Override
    public String getToolTipText() {
        return "Add reference allele to HDF5 genotypes";
    }
}
