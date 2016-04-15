/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import java.util.ArrayList;
import java.util.Arrays;

import com.google.common.primitives.Ints;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 *
 * @author qs24
 */
public class VCFUtil {
    // variables for calculating OS and PL for VCF, might not be in the correct class
    private static double error;
    private static double v1;
    private static double v2;
    private static double v3;
    private static int[][][] myGenoScoreMap;
    public static final int VCF_DEFAULT_MAX_NUM_ALLELES = 3;

    static  {
        error = 0.001; //TODO this seems low, is this the standard
        v1 = Math.log10(1.0 - error * 3.0 /4.0);
        v2 = Math.log10(error/4);
        v3 = Math.log10(0.5 - (error/4.0));
        myGenoScoreMap = new int[128][128][];
        for (int i = 0; i < 128; i++) {
            for (int j = 0; j < 128; j++) {
                myGenoScoreMap[i][j]= calcScore(i, j);
            }
        }
    }

    private VCFUtil ()
    {
        
    }
    

    
    public static int[] getScore(int i, int j) {
        if(i>127 || j>127) return calcScore(i,j);
        return myGenoScoreMap[i][j];
    }
    
    // Calculate QS and PL for VCF might not be in the correct class
    private static int[] calcScore (int a, int b)
    {   
        int[] results= new int[4];
        int n = a + b;
        int m = a;
        if (b > m) {
            m = b;
        }

        double fact = 0;
        if (n > m) {
            for (int i = n; i > m; i--) {
               fact += Math.log10(i);
            }
            for (int i = 1; i <= (n - m); i++){
               fact -= Math.log10(i);
            }
        }
        double aad = Math.pow(10, fact + (double)a * v1 + (double)b * v2);
        double abd = Math.pow(10, fact + (double)n * v3);
        double bbd = Math.pow(10, fact + (double)b * v1 + (double)a * v2);
        double md = aad;
        if (md < abd) {
            md = abd;
        }
        if (md < bbd) {
            md = bbd;
        }
        int gq = 0;
        if ((aad + abd + bbd) > 0) {
            gq = (int)(md / (aad + abd + bbd) * 100);
        }
        
        int aa =(int) (-10 * (fact + (double)a * v1 + (double)b * v2));
        int ab =(int) (-10 * (fact + (double)n * v3));
        int bb =(int) (-10 * (fact + (double)b * v1 + (double)a * v2));
        
        m = aa;
        if (m > ab) {
            m = ab;
        }
        if (m>bb) {
            m = bb;
        }
        aa -= m;
        ab -= m;
        bb -= m;
        results[0] = aa > 255 ? 255 : aa;
        results[1] = ab > 255 ? 255 : ab;
        results[2] = bb > 255 ? 255 : bb;
        results[3] = gq;
        
        return results;
    }
    
     public static byte resolveVCFGeno(byte[] alleles, int[][] allelesInTaxa, int tx) {
        int[] alleleDepth = new int[allelesInTaxa.length];
        for (int i=0; i<allelesInTaxa.length; i++)
        {
            alleleDepth[i] = allelesInTaxa[i][tx];
        }
        return resolveVCFGeno(alleles, alleleDepth);
    }
     
     public static byte resolveVCFGeno(byte[] alleles, int[] alleleDepth) { 
        int depth = 0;
        for (int i = 0; i < alleleDepth.length; i++) {
            depth += alleleDepth[i];
        }
        if (depth == 0) {
            return (byte)((GenotypeTable.UNKNOWN_ALLELE << 4) | GenotypeTable.UNKNOWN_ALLELE);
        }
        int max = 0;
        byte maxAllele = GenotypeTable.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = GenotypeTable.UNKNOWN_ALLELE;
        for (int i = 0; i < alleles.length; i++) {
            if (alleleDepth[i] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = alleleDepth[i];
                maxAllele = alleles[i];
            } else if (alleleDepth[i] > nextMax) {
                nextMax = alleleDepth[i];
                nextMaxAllele = alleles[i];
            }
        }
        if (alleles.length == 1) {
            return (byte)((alleles[0] << 4) | alleles[0]);
        } else {
            max = (max > 127) ? 127 : max;
            nextMax = (nextMax > 127) ? 127 : nextMax;
            int[] scores = getScore(max, nextMax);
            if ((scores[1] <= scores[0]) && (scores[1] <= scores[2])) {
                return (byte)((maxAllele << 4) | nextMaxAllele);
            } else if ((scores[0] <= scores[1]) && (scores[0] <= scores[2])) {
                return (byte)((maxAllele << 4) | maxAllele);
            } else {
                return (byte)((nextMaxAllele << 4) | nextMaxAllele);
            }
        }
     }
     

     public static int[] resolveRefSorted(int[] sortedAlleles, byte refAllele) {
         int[] sortedAllelesResolved = new int[sortedAlleles.length];
         int indexOfRefAllele = Ints.indexOf(sortedAlleles, refAllele);
         
         //If indexOfRefAllele is -1 the refAllele is not contained in sortedAlleles
         //We need to add it
         if (indexOfRefAllele < 0) {
             //Set index to 0 as we want the Ref in the first position
             indexOfRefAllele = 0;
             
             if(refAllele != GenotypeTable.UNKNOWN_ALLELE) {
                 int[] sortedAllelesExpanded = new int[sortedAlleles.length+1];
                 sortedAllelesExpanded[0] = refAllele;
                 for(int i = 0; i<sortedAlleles.length; i++) {
                     sortedAllelesExpanded[i+1] = sortedAlleles[i];
                 }
                 sortedAllelesResolved = sortedAllelesExpanded;
             }
             else {
                 for(int i = 0; i < sortedAllelesResolved.length; i++) {
                     sortedAllelesResolved[i] = sortedAlleles[i];
                 }
             }
         }
         //Resort sorted Alleles if ref is not first in array
         else if (indexOfRefAllele != 0) {
             sortedAllelesResolved[0] = refAllele;
             
             for(int i = indexOfRefAllele; i>0; i--) { 
                 sortedAllelesResolved[i] = sortedAlleles[i-1]; 
             }
             for(int i = indexOfRefAllele+1; i<sortedAllelesResolved.length;i++) {
                 sortedAllelesResolved[i] = sortedAlleles[i];
             }
         }
         else {
             for(int i = 0; i < sortedAllelesResolved.length;i++) {
                 sortedAllelesResolved[i] = sortedAlleles[i];
             }
         }
         
         return sortedAllelesResolved;
     }
     
     public static boolean indelInKnownVariant(String[] knownVariants) {
         for(String variant:knownVariants) {
             if(variant.length()>1) {
                 return true;
             }
         }
         return false;
     }
     
     public static Tuple<int[], String[]> resolveSortedAndKnownVariantsExport(int[] sortedAllelesInput, String[] knownVariantsInput) {
         int[] sortedAlleles = Arrays.copyOf(sortedAllelesInput, sortedAllelesInput.length);
         String[] knownVariants = Arrays.copyOf(knownVariantsInput,knownVariantsInput.length);
         if(knownVariants.length>0) {
             //ReOrder based on variant alleles
             //Store a tempSortedAlleles so we can appropriately handle hapmap to vcf
             //int[] tempSortedAlleles = new int[knownVariants.length];
             
             //ArrayList to hold the Sorted Alleles Indices Temporarily as the ordering will change
             ArrayList<Integer> tempSortedAlleles = new ArrayList<Integer>();
             
             //Loop through all the knownVariants and check to see if we have an indel
             boolean knownVariantIndel = VCFUtil.indelInKnownVariant(knownVariants);
             
             //If we do have an indel, we can add the variants after picking off the first character to the tempSortedAlleles
             if(knownVariantIndel) {
                 //Loop through the variants
                 for(int i = 0; i < knownVariants.length; i++) {
                     //Pull off the first character if it exists
                     if(knownVariants[i].length()>1) {
                         String parsedVariant = knownVariants[i].substring(1);
                         tempSortedAlleles.add((int)NucleotideAlignmentConstants.getNucleotideAlleleByte(parsedVariant.charAt(0)));
                     }
                     else {
                         //Mark as deletion
                         tempSortedAlleles.add((int)NucleotideAlignmentConstants.getNucleotideAlleleByte('-'));
                     }
                 }
             } else {
                 //If we dont have an indel, we can add it to the allele array
                 if(sortedAlleles.length<knownVariants.length){
                     //Clear it out, we probably dont need to do this
                     tempSortedAlleles = new ArrayList<Integer>();
                 }
                 int nIndex = -1;
                 for(int i = 0; i<knownVariants.length; i++) {
                     //ZRM22 Mar 22
                     if(knownVariants[i].charAt(0)!='N') {
                         tempSortedAlleles.add((int)NucleotideAlignmentConstants.getNucleotideAlleleByte(knownVariants[i].charAt(0)));
                     }
                     else {
                         //If N is in our known Variants list but we do not have an indel, we need to remove it
                         nIndex = i;
                     }
                 }
                 if(nIndex != -1) {
                     //if we have an N we need to resize KnownVariants
                     String[] knownVariantsSmall = new String[knownVariants.length-1];
                     for(int i = 0; i<knownVariants.length; i++) {
                         if(i < nIndex) {
                             knownVariantsSmall[i] = knownVariants[i];
                         }
                         else if(i > nIndex) {
                             knownVariantsSmall[i-1] = knownVariants[i];
                         }
                     }
                     knownVariants = knownVariantsSmall;
                 }
             }
             //END ZRM22 Jan7
             
             //Make a copy of KnownVaraints in case we need to add some
             ArrayList<String> knownVariantsList = new ArrayList<String>();
             boolean indelsExist = false;
             boolean indelsInKnownVariants = VCFUtil.indelInKnownVariant(knownVariants);
             if(indelsInKnownVariants) {
                 indelsExist = true;
             }
             
             //Go through sorted alleles and also check for indels
             for(int i = 0 ;i<sortedAlleles.length; i++) {
                 if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]).equals("-")) {
                     indelsExist = true;
                 }
             }
             //Move To Function/
             
             for(String variant:knownVariants) {
                 if(indelsExist && !indelsInKnownVariants) {
                     knownVariantsList.add("N"+variant);
                 }
                 else {
                     knownVariantsList.add(variant);
                 }
             }
             
             //Go through sorted alleles
             for(int i = 0 ;i<sortedAlleles.length; i++) {
             //If a sorted allele is not in tempSortedAlleles,
                 if(!tempSortedAlleles.contains(sortedAlleles[i])) {
                     //if its not add it to sorted alleles and knownVariants
                     tempSortedAlleles.add(sortedAlleles[i]);
                     //Check for an indel
                     if(indelsExist) {
                         if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]).equals("-")) {
                             knownVariantsList.add("N");
                         }
                         else {
                             knownVariantsList.add("N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
                         }
    //                     knownVariantsList.add("N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
                     }
                     else {
                         knownVariantsList.add(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
                     }
                 }
             }
             //reset knownVariants and sortedAlleles to reflect the changes
             String[] knownVariantsExtended = new String[knownVariantsList.size()];
             for(int i = 0; i < knownVariantsExtended.length; i++) {
                 knownVariantsExtended[i] = knownVariantsList.get(i);
             }
             knownVariants = knownVariantsExtended;
             
             int[] sortedAllelesExtended = new int[tempSortedAlleles.size()];
             for(int i = 0; i < sortedAllelesExtended.length; i++) {
                 sortedAllelesExtended[i] = tempSortedAlleles.get(i);
             }
             sortedAlleles = sortedAllelesExtended;
             //sortedAlleles = tempSortedAlleles.toArray(new int[tempSortedAlleles.size()]);
         }
         else {
             //No known variants, but we need to handle indels
             int indelIndex = -1;
             //loop through sorted alleles
             for(int i = 0; i<sortedAlleles.length; i++) {
                 //if we find an indel mark the index and set a boolean
                 if(sortedAlleles[i] == (int)NucleotideAlignmentConstants.getNucleotideAlleleByte('-')) {
                     indelIndex = i;
                     break;
                 }
             }
             
             knownVariants = new String[sortedAlleles.length];
             for(int i = 0; i<knownVariants.length; i++) {
                 if(indelIndex==-1) {
                     knownVariants[i] = ""+ NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]);
                 }
                 else {
                     if(indelIndex == i) {
                         knownVariants[i] = "N";
                     }
                     else {
                         knownVariants[i] = "N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]);
                     }
                 }
             }
         }
         return new Tuple<int[],String[]>(sortedAlleles,knownVariants);         
     }
}
