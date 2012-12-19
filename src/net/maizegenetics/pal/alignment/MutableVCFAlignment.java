/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.alignment;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author yz79
 */
public class MutableVCFAlignment extends MutableNucleotideAlignment implements MutableAlignment {
    
    // allelic depth information, byte[allele][taxa][site] = depth, max: 255
    private byte[][][] myAlleleDepth;
    
    // possible alleles for each site
    private byte[][] myCommonAlleles;
    
    private final int myMaxNumAlleles = 3;
    
    private MutableVCFAlignment(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a, maxNumTaxa, maxNumSites);
        initAllelicDepthArrays(maxNumTaxa, maxNumSites);
    }
    
    private MutableVCFAlignment(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        super(idGroup, initNumSites, maxNumTaxa, maxNumSites);
        initAllelicDepthArrays(maxNumTaxa, maxNumSites);
    }
    
    private MutableVCFAlignment(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        super(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
        initAllelicDepthArrays(idGroup.size(), siteNames.length);
    }
    
    public static MutableVCFAlignment getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        if (a.getAlleleEncodings() == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES) {
            return new MutableVCFAlignment(a, maxTaxa, maxNumSites);
        } else {
            throw new IllegalArgumentException("MutableNucleotideAlignment: getInstance: alignment must be nucleotide data.");
        }
    }
    
    public static MutableVCFAlignment getInstance(IdGroup idGroup, int maxNumSites) {
        return new MutableVCFAlignment(idGroup, 0, idGroup.getIdCount(), maxNumSites);
    }
    
    public static MutableVCFAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        return new MutableVCFAlignment(idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }
    
    public static MutableVCFAlignment getInstance(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        return new MutableVCFAlignment(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
    }
    
    public void initAllelicDepthArrays(int numMaxTaxa, int numMaxSites) {
        myAlleleDepth = new byte[myMaxNumAlleles][numMaxTaxa][numMaxSites];
        myCommonAlleles = new byte[myMaxNumAlleles][numMaxSites];
        for (int i = 0; i < myAlleleDepth.length; i++) {
            for (int j = 0; j < myAlleleDepth[i].length; j++) {
                for (int k = 0; k < myAlleleDepth[i][j].length; k++) {
                    myAlleleDepth[i][j][k] = (byte) -1;
                }
                myCommonAlleles[i][j] = (byte) -1;
            }
        }
    }
    
    @Override
    public byte[] getAlleles(int site) {
        ArrayList<Byte> outArray = new ArrayList<Byte>();
        for (int i = 0; i < myCommonAlleles.length; i++) {
            if (myCommonAlleles[i][site] != (byte) -1) {
                outArray.add(myCommonAlleles[i][site]);
            } else {
                break;
            }
        }
        byte[] out = new byte[outArray.size()];
        for (int i = 0; i < outArray.size(); i++) {
            out[i] = outArray.get(i).byteValue();
        }
        return out;
    }
    
    @Override
    public byte[] getDepthForAllele(int taxon, int site) {
        ArrayList<Byte> outArray = new ArrayList<Byte>();
        for (int i = 0; i < myAlleleDepth.length; i++) {
            if (myAlleleDepth[i][taxon][site] != (byte) -1) {
                outArray.add(myAlleleDepth[i][taxon][site]);
            } else {
                break;
            }
        }
        byte[] out = new byte[outArray.size()];
        for (int i = 0; i < outArray.size(); i++) {
            out[i] = outArray.get(i).byteValue();
        }
        return out;
    }
    
    public void setCommonAllele(int site, byte[] values) {
        for (int i = 0; i < values.length; i++) {
            myCommonAlleles[i][site] = values[i];
        }
    }
    
    public void setDepthForAllele(int taxon, int site, byte[] values) {
        for (int i = 0; i < values.length; i++) {
            myAlleleDepth[i][taxon][site] = values[i];
        }
    }
}
