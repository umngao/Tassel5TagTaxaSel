/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import net.maizegenetics.gwas.imputation.EmissionProbability;
import net.maizegenetics.gwas.imputation.TransitionProbability;
import net.maizegenetics.gwas.imputation.ViterbiAlgorithm;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;


/**
 * Imputation methods that relies on an alignment of possible gametes (donorAlign), and optimized for
 * highly heterozygous samples.  It works at the scale of 64 sites to accelerate searching
 * through all possible parental combinations for a window.  
 * 
 * It only begins each 64 site block and expand outward to find set number of minor
 * alleles in the target sequence, and then it looks for all possible parents.
 * 
 * TODO: remove bad sites, add hets back in
 * 4. Viterbi between sub-chromosome segments
 * 5. Pre-compute IBD probabilities for different number of sites and divergences 
 * (Doing this by poisson on the fly is too slow, and use this rank best hybrid or inbred
 * 6. Error rates by taxon
 * 
 * TESTS:
 * 1.  Do not impute the invariant sites - tested doesn't matter
 * 2.  Change imputation to fixed windows - minor allele is a little more sensitive, but not by much
 * 3.  Change resolve regionally to HMM - it may work better with 2
 * 
 * @author edbuckler, kswarts, aromero
 */
public class MinorWindowViterbiImputation {
    private Alignment unimpAlign;  //the unimputed alignment to be imputed, unphased
    private Alignment donorAlign;  //these are the reference haplotypes, that must be homozygous
    private int testing=0;  //level of reporting to stdout
    //major and minor alleles can be differ between the donor and unimp alignment 
    //these Bit sets keep track of the differences
    private OpenBitSet swapMjMnMask, invariantMask, errorMask, goodMask;
    private boolean isSwapMajorMinor=false;  //if swapped try to fix it
    
    //variables for tracking accuracy 
    private int totalRight=0, totalWrong=0, totalHets=0; //global variables tracking errors on the fly
    private int[] siteErrors, siteCallCnt, taxonErrors, taxonCallCnt;  //error recorded by sites
    
    
    private int blocks=-1;  //total number 64 site words in the alignment
    
    private int minMajorRatioToMinorCnt=10;  //refinement of minMinorCnt to account for regions with all major
    private int maxDonorHypotheses=10;  //number of hypotheses of record from an inbred or hybrid search of a focus block
    
    private double maximumInbredError=0.02;  //inbreds are tested first, if too much error hybrids are tested.
  //  private double maxViterbiErro
    private int minTestSites=100;  //minimum number of compared sites to find a hit
    
    //matrix to hold divergence comparisons between the target taxon and donor haplotypes for each block
    //dimensions: [donor taxa][sites, same count, diff count, het count][blocks] 
    private byte[][][] allDist;
    
    //initialize the transition matrix
    double[][] transition = new double[][] {
                    {.999,.0001,.0003,.0001,.0005},
                    {.0002,.999,.00005,.00005,.0002},
                    {.0002,.00005,.999,.00005,.0002},
                    {.0002,.00005,.00005,.999,.0002},
                    {.0005,.0001,.0003,.0001,.999}
        };
    //initialize the emission matrix, states (5) in rows, observations (3) in columns
    double[][] emission = new double[][] {
                    {.98,.001,.001},
                    {.6,.2,.2},
                    {.4,.2,.4},
                    {.2,.2,.6},
                    {.001,.001,.98}
    };
    TransitionProbability tp = new TransitionProbability();
    EmissionProbability ep = new EmissionProbability();
    
    int[] donorIndices;
       
    private static int highMask=0xF0;
    private static int lowMask=0x0F;
    
    

    /**
     * 
     * @param donorFile should be phased haplotypes 
     * @param unImpTargetFile sites must match exactly with donor file
     * @param exportFile output file of imputed sites
     * @param minMinorCnt determines the size of the search window, low recombination 20-30, high recombination 10-15
     * @param hybridMode true=do hybrid search
     */
    public MinorWindowViterbiImputation(String donorFile, String unImpTargetFile, 
            String exportFile, int minMinorCnt, int minTestSites, int minSitesPresent, 
            double maxHybridErrorRate, boolean hybridMode, boolean imputeDonorFile) {
        this.minTestSites=minTestSites;
        donorAlign=ImportUtils.readFromHapmap(donorFile, false, (ProgressListener)null);
        donorAlign.optimizeForTaxa(null);
        System.out.printf("Donor taxa:%d sites:%d %n",donorAlign.getSequenceCount(),donorAlign.getSiteCount());        
        donorIndices=new int[donorAlign.getSequenceCount()];
        for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i;}
        if(imputeDonorFile) {
            unimpAlign=donorAlign;
        } else {
            unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
            unimpAlign.optimizeForTaxa(null);
        }
        siteErrors=new int[unimpAlign.getSiteCount()];
        siteCallCnt=new int[unimpAlign.getSiteCount()];
        taxonErrors=new int[unimpAlign.getSequenceCount()];
        taxonCallCnt=new int[unimpAlign.getSequenceCount()];
        tp.setTransitionProbability(transition);
        ep.setEmissionProbability(emission);
        System.out.printf("Unimputed taxa:%d sites:%d %n",unimpAlign.getSequenceCount(),unimpAlign.getSiteCount());
        
        System.out.println("Creating mutable alignment");
        blocks=unimpAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        createMaskForAlignmentConflicts();
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.unimpAlign);  //output data matrix
        long time=System.currentTimeMillis();
        int countFullLength=0;
        for (int taxon = 0; taxon < unimpAlign.getSequenceCount(); taxon+=1) {
            int blocksSolved=0;
            if(imputeDonorFile){
                donorIndices=new int[donorAlign.getSequenceCount()-1];
                for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i; if(i>=taxon) donorIndices[i]++;}
            }
            DonorHypoth[][] regionHypth=new DonorHypoth[blocks][maxDonorHypotheses];
            String name=unimpAlign.getIdGroup().getIdentifier(taxon).getFullName();
            BitSet[] maskedTargetBits=arrangeMajorMinorBtwAlignments(unimpAlign,taxon);
            System.out.printf("Imputing %d:%s Mj:%d, Mn:%d Unk:%d ... ", taxon,name,maskedTargetBits[0].cardinality(),
                    maskedTargetBits[1].cardinality(), countUnknown(mna,taxon));
            if(unimpAlign.getTotalNotMissingForTaxon(taxon)<minSitesPresent) continue;
            calcInbredDist(maskedTargetBits);
            for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
                if(resultRange==null) continue; //no data in the focus Block
                //search for the best inbred donors for a segment
                regionHypth[focusBlock]=getBestInbredDonors(taxon, resultRange[0],resultRange[2], focusBlock, donorIndices);
            }
            boolean foundit=apply1or2Haplotypes(taxon, regionHypth,  mna, maskedTargetBits, maxHybridErrorRate);
            if(foundit) {
                countFullLength++;
            } else {
                blocksSolved=solveRegionally3(mna, taxon, regionHypth, hybridMode, maskedTargetBits, 
                        minMinorCnt, maxHybridErrorRate);
            }
            int unk=countUnknown(mna,taxon);
            System.out.printf("Viterbi: %s BlocksSolved:%d Done Unk: %d PropMissing:%g ", foundit, blocksSolved, unk, (double)unk/(double)mna.getSiteCount());
            double errRate=(double)taxonErrors[taxon]/(double)(taxonCallCnt[taxon]+taxonErrors[taxon]); 
            System.out.printf("ErR:%g %n", errRate);
        }
        System.out.println("");
        System.out.println("Time:"+(time-System.currentTimeMillis()));
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        System.out.println(s.toString());
        double errRate=(double)totalWrong/(double)(totalRight+totalWrong);
        System.out.printf("TotalRight %d  TotalWrong %d TotalHets: %d ErrRateExcHet:%g %n",totalRight, totalWrong, totalHets, errRate);
//        for (int i = 0; i < siteErrors.length; i++) {
//            System.out.printf("%d %d %d %g %g %n",i,siteCallCnt[i],siteErrors[i], 
//                    (double)siteErrors[i]/(double)siteCallCnt[i], unimpAlign.getMinorAlleleFrequency(i));
//        }
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
        System.out.printf("%d %g %d %n",minMinorCnt, maximumInbredError, maxDonorHypotheses);
        System.out.printf("Full length imputed %d %n ",countFullLength);
        
    }

    private boolean apply1or2Haplotypes(int taxon, DonorHypoth[][] regionHypth, MutableNucleotideAlignment mna,
            BitSet[] maskedTargetBits, double maxHybridErrorRate) {
        
        //do flanking search 
        if(testing==1) System.out.println("Starting complete hybrid search");
        int[] d=getAllBestDonorsAcrossChromosome(regionHypth,5);
        DonorHypoth[] best2donors=getBestHybridDonors(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, d, d, true);
        if(testing==1) System.out.println(Arrays.toString(best2donors));
        ArrayList<DonorHypoth> goodDH=new ArrayList<DonorHypoth>();
        for (DonorHypoth dh : best2donors) {
            if((dh!=null)&&(dh.getErrorRate()<maxHybridErrorRate)) {
                if(dh.isInbred()==false){
                    dh=getStateBasedOnViterbi(dh);
                }
                if(dh!=null) goodDH.add(dh);
            }
        }
        if(goodDH.size()==0) return false;
        DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
        for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
        setAlignmentWithDonors(vdh,false,mna);
        return true;
    }    
    
    
    /**
     * Create mask for all sites where major & minor are swapped in alignments
     */
    private void createMaskForAlignmentConflicts() {
        goodMask=new OpenBitSet(unimpAlign.getSiteCount());
        errorMask=new OpenBitSet(unimpAlign.getSiteCount());
        swapMjMnMask=new OpenBitSet(unimpAlign.getSiteCount());
        invariantMask=new OpenBitSet(unimpAlign.getSiteCount());
        int siteConflicts=0, swaps=0, invariant=0, good=0;
        for (int i = 0; i < unimpAlign.getSiteCount(); i++) {
            /*we have three classes of data:  invariant in one alignment, conflicts about minor and minor,
            *swaps of major and minor.  Adding the invariant reduces imputation accuracy.
            *the major/minor swaps should be flipped in the comparisons
            */
            if(donorAlign.getMinorAllele(i)==Alignment.UNKNOWN_ALLELE) {
                invariant++;
                invariantMask.set(i);
                goodMask.set(i);
            } else 
            if((donorAlign.getMajorAllele(i)==unimpAlign.getMinorAllele(i))&&(donorAlign.getMinorAllele(i)==unimpAlign.getMajorAllele(i))) {
                swaps++;
                swapMjMnMask.set(i);
                goodMask.set(i);
            } else
            if((donorAlign.getMajorAllele(i)!=unimpAlign.getMajorAllele(i))) {
                siteConflicts++;
                errorMask.set(i);
                goodMask.set(i);
            } 
            
        }
        goodMask.not();
        System.out.println("invariant in donor:"+invariant+" swapConflicts:"+swaps+" errors:"+siteConflicts);
    }
    
    private int countUnknown(Alignment a, int taxon) {
        int cnt=0;
        for (int i = 0; i < a.getSiteCount(); i++) {
            if(a.getBase(taxon, i)==Alignment.UNKNOWN_DIPLOID_ALLELE) cnt++;
        }
        return cnt;
    }

    
    private BitSet[] arrangeMajorMinorBtwAlignments(Alignment unimpAlign, int bt) { 
        OpenBitSet mjTbs=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 0));
        OpenBitSet mnTbs=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 1));
        mjTbs.and(goodMask);
        mnTbs.and(goodMask);
        if(isSwapMajorMinor) {
            OpenBitSet newmj=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 1));
            OpenBitSet newmn=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 0));
            newmj.and(swapMjMnMask);
            newmn.and(swapMjMnMask);
            mjTbs.or(newmj);
            mnTbs.or(newmn);
        }
        BitSet[] result={mjTbs,mnTbs};
        return result;
    }
    
    
    /**
     * If the target regions has Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * TODO:  This can be done much more robustly.
     * @param mna
     * @param targetTaxon
     * @param regionHypth 
     */
    private void solveRegionally2(MutableNucleotideAlignment mna, int targetTaxon, 
            DonorHypoth[][] regionHypth, boolean hybridMode, BitSet[] maskedTargetBits, int minMinorCnt) {
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
            if(resultRange==null) continue; //no data in the focus Block
             //search for the best hybrid donors for a segment
            if(hybridMode&&(regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()>maximumInbredError)) {
                long[] mjTRange=maskedTargetBits[0].getBits(resultRange[0],resultRange[2]);
                long[] mnTRange=maskedTargetBits[1].getBits(resultRange[0],resultRange[2]);
                regionHypth[focusBlock]=getBestHybridDonors(targetTaxon, mjTRange, mnTRange, resultRange[0],resultRange[2], 
                        focusBlock, regionHypth[focusBlock], false);
                //hybrid++;
            }
        }
        System.out.println("");
        System.out.print("Donors:");
        DonorHypoth[][] theLH=new DonorHypoth[blocks][];
        theLH[0]=regionHypth[0];
        DonorHypoth[][] theRH=new DonorHypoth[blocks][];
        theRH[blocks-1]=regionHypth[blocks-1];
        for (int focusBlock = 1; focusBlock < blocks; focusBlock++) {
            if((regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()<this.maximumInbredError)) {
                    theLH[focusBlock]=regionHypth[focusBlock];   
                } else {theLH[focusBlock]=theLH[focusBlock-1];}
        }
        for (int focusBlock = blocks-2; focusBlock >=0; focusBlock--) {
            if((regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()<this.maximumInbredError)) {
                    theRH[focusBlock]=regionHypth[focusBlock];   
                } else {theRH[focusBlock]=theRH[focusBlock+1];}
        }
        int focusBlockdone=0;
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            if(theLH[focusBlock][0]==null) continue;
            if(theLH[focusBlock]==theRH[focusBlock]) {
//            if(theLH[focusBlock][0]!=null) {
                DonorHypoth[] newHybrid=new DonorHypoth[1];
                newHybrid[0]=new DonorHypoth(theLH[focusBlock][0].targetTaxon,theLH[focusBlock][0].donor1Taxon, 
                        theLH[focusBlock][0].donor1Taxon, theLH[focusBlock][0].startBlock, focusBlock, theLH[focusBlock][0].endBlock);
                setAlignmentWithDonors(newHybrid,true,mna);
                focusBlockdone++;
                //System.out.print(" "+focusBlock+":"+theLH[focusBlock][0].donor1Taxon+":"+theLH[focusBlock][0].donor2Taxon);

            } 
            else {
                if(theRH[focusBlock][0]==null) continue;
                if(theLH[focusBlock][0].donor1Taxon>=-1) continue;
                DonorHypoth[] newHybrid=new DonorHypoth[1];
                newHybrid[0]=new DonorHypoth(theLH[focusBlock][0].targetTaxon,theLH[focusBlock][0].donor1Taxon, 
                        theRH[focusBlock][0].donor1Taxon, theLH[focusBlock][0].startBlock, focusBlock, theRH[focusBlock][0].endBlock);
                setAlignmentWithDonors(newHybrid,true,mna);
                focusBlockdone++;
               //System.out.print(" "+focusBlock+":"+newHybrid[0].donor1Taxon+":"+newHybrid[0].donor2Taxon);
            }

        }
        int leftNullCnt=0;
        for (DonorHypoth[] dh : theLH) {if(dh[0]==null) leftNullCnt++;}
        if(testing==1) System.out.printf("block:%d focusBlockdone:%d nullLeft:%d %n",blocks,focusBlockdone, leftNullCnt);
    }
    
    /**
     * If the target regions has Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * TODO:  This can be done much more robustly.
     * @param mna
     * @param targetTaxon
     * @param regionHypth 
     */ 
    private int solveRegionally3(MutableNucleotideAlignment mna, int targetTaxon, 
            DonorHypoth[][] regionHypth, boolean hybridMode, BitSet[] maskedTargetBits, 
            int minMinorCnt, double maxHybridErrorRate) {
        int focusBlockdone=0, inbred=0, hybrid=0, noData=0;
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            if((regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()<maximumInbredError)) {
                setAlignmentWithDonors(regionHypth[focusBlock],true,mna);
                focusBlockdone++; 
                inbred++;
            }  else {
                int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),
                        maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
                if(resultRange==null) {noData++; continue; }//no data in the focus Block
                 //search for the best hybrid donors for a segment
                if(hybridMode&&(regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()>maximumInbredError)) {
                    long[] mjTRange=maskedTargetBits[0].getBits(resultRange[0],resultRange[2]);
                    long[] mnTRange=maskedTargetBits[1].getBits(resultRange[0],resultRange[2]);
                    regionHypth[focusBlock]=getBestHybridDonors(targetTaxon, mjTRange, mnTRange, resultRange[0],resultRange[2], 
                            focusBlock, regionHypth[focusBlock], false);
                }
          //      if((regionHypth[focusBlock][0]!=null)) System.out.printf("%d %d %g %n",targetTaxon, focusBlock, regionHypth[focusBlock][0].getErrorRate() );
                if((regionHypth[focusBlock][0]!=null)&&(regionHypth[focusBlock][0].getErrorRate()<maxHybridErrorRate)) {
                    setAlignmentWithDonors(regionHypth[focusBlock],true,mna);
                    focusBlockdone++;
                    hybrid++;
                } else {
                    if(regionHypth[focusBlock][0]==null) {noData++;} 
//                    else{
//                        System.out.println(regionHypth[focusBlock][0].getErrorRate());
//                    }
                }
            }
        }
        
        int leftNullCnt=0;
        for (DonorHypoth[] dh : regionHypth) {if(dh[0]==null) leftNullCnt++;}
        if(testing==1) System.out.printf("targetTaxon:%d hybridError:%g block:%d focusBlockdone:%d null:%d inbredDone:%d hybridDone:%d noData:%d %n",
                targetTaxon, maxHybridErrorRate, blocks,focusBlockdone, leftNullCnt, inbred, hybrid, noData);
        return focusBlockdone;
    }
    
    private DonorHypoth getStateBasedOnViterbi(DonorHypoth dh) {
        int startSite=dh.startBlock*64;
        int endSite=(dh.endBlock*64)+63;
        if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
        int sites=endSite-startSite+1;
        byte[] calls=new byte[sites];
        //System.out.printf("%d %d %d %n",dh.donor1Taxon, startSite, endSite+1);
        byte[] d1b=donorAlign.getBaseRange(dh.donor1Taxon, startSite, endSite+1);
        byte[] d2b=donorAlign.getBaseRange(dh.donor2Taxon, startSite, endSite+1);
        byte[] t1b=unimpAlign.getBaseRange(dh.targetTaxon, startSite, endSite+1);
        int informSites=0, nonMendel=0;
        ArrayList<Byte> nonMissingObs = new ArrayList<Byte>();
        ArrayList<Integer> snpPositions = new ArrayList<Integer>();
        for(int cs=0; cs<sites; cs++) {
            if(t1b[cs]==Alignment.UNKNOWN_DIPLOID_ALLELE) continue; 
            if(d1b[cs]==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
            if(d2b[cs]==Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
            if(t1b[cs]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) continue; 
            if(d1b[cs]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) continue;
            if(d2b[cs]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) continue;
            if(d1b[cs]==d2b[cs]) {
                if(t1b[cs]!=d1b[cs]) nonMendel++;
                continue;
            }
            informSites++;
            byte state=1;
            if(t1b[cs]==d1b[cs]) {state=0;}
            else if(t1b[cs]==d2b[cs]) {state=2;}
            nonMissingObs.add(state);
            snpPositions.add(cs+startSite);
        }
        if(informSites<10) return null;
        double nonMendRate=(double)nonMendel/(double)informSites;
        if(testing==1) System.out.printf("NonMendel:%d InformSites:%d ErrorRate:%g %n",nonMendel, informSites, nonMendRate);
        if(nonMendRate>this.maximumInbredError*5) return null;
        byte[] informStates=new byte[informSites];
        for (int i = 0; i < informStates.length; i++) informStates[i]=nonMissingObs.get(i);
        int[] pos=new int[informSites];
        for (int i = 0; i < pos.length; i++) pos[i]=snpPositions.get(i);
        int chrlength = donorAlign.getPositionInLocus(endSite) - donorAlign.getPositionInLocus(startSite);
        tp.setAverageSegmentLength( chrlength / sites );
        tp.setPositions(pos);
        
	double probHeterozygous=0.5;
        double phom = (1 - probHeterozygous) / 2;
        double[] pTrue = new double[]{phom, .25*probHeterozygous ,.5 * probHeterozygous, .25*probHeterozygous, phom};
		
        ViterbiAlgorithm va = new ViterbiAlgorithm(informStates, tp, ep, pTrue);
	va.calculate();
        if(testing==1) System.out.println("Input:"+Arrays.toString(informStates));
        byte[] resultStates=va.getMostProbableStateSequence();
        if(testing==1) System.out.println("Resul:"+Arrays.toString(resultStates));
        DonorHypoth dh2=new DonorHypoth(dh.targetTaxon,dh.donor1Taxon, 
                        dh.donor2Taxon, dh.startBlock, dh.focusBlock, dh.endBlock);
        int currPos=0;
        for(int cs=0; cs<sites; cs++) {
            calls[cs]=(resultStates[currPos]==1)?(byte)1:(byte)(resultStates[currPos]/2); //converts the scale back to 0,1,2 from 0..4
            if((pos[currPos]<cs+startSite)&&(currPos<resultStates.length-1)) currPos++;
        }
        dh2.phasedResults=calls;
        return dh2;
    }
    
    
    
    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param focusBlock
     * @param minMinorCnt
     * @return arrays of blocks {startBlock, focusBlock, endBlock}
     */
    private int[] getBlockWithMinMinorCount(long[] mjT, long[] mnT, int focusBlock, int minMinorCnt) {
        int majorCnt=Long.bitCount(mjT[focusBlock]);
        int minorCnt=Long.bitCount(mnT[focusBlock]);
        int endBlock=focusBlock, startBlock=focusBlock;
        int minMajorCnt=minMinorCnt*minMajorRatioToMinorCnt;
        while((minorCnt<minMinorCnt)&&(majorCnt<minMajorCnt)) {
            boolean preferMoveStart=(focusBlock-startBlock<endBlock-focusBlock)?true:false;
            if(startBlock==0) {preferMoveStart=false;}
            if(endBlock==blocks-1) {preferMoveStart=true;}
            if((startBlock==0)&&(endBlock==blocks-1)) break;
            if(preferMoveStart) {//expand start
                startBlock--;
                minorCnt+=Long.bitCount(mnT[startBlock]);
                majorCnt+=Long.bitCount(mjT[startBlock]);
            } else { //expand end
                endBlock++;
                minorCnt+=Long.bitCount(mnT[endBlock]); 
                majorCnt+=Long.bitCount(mjT[startBlock]);
            } 
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
    }
  
    
    /**
     * Calculates an inbred distance
     * TODO:  Should be converted to the regular same, diff, hets, & sum
     * @param modBitsOfTarget
     * @return array with sites, same count, diff count, het count
     */
    private void calcInbredDist(BitSet[] modBitsOfTarget) {
        allDist=new byte[donorAlign.getSequenceCount()][4][blocks];
        long[] iMj=modBitsOfTarget[0].getBits();
        long[] iMn=modBitsOfTarget[1].getBits();
        for (int donor1 = 0; donor1 < allDist.length; donor1++) {
            long[] jMj=donorAlign.getAllelePresenceForAllSites(donor1, 0).getBits();
            long[] jMn=donorAlign.getAllelePresenceForAllSites(donor1, 1).getBits();
            for (int i = 0; i <blocks; i++) {
                long same = (iMj[i] & jMj[i]) | (iMn[i] & jMn[i]);
                long diff = (iMj[i] & jMn[i]) | (iMn[i] & jMj[i]);
                long hets = same & diff;
                int sameCnt = BitUtil.pop(same);
                int diffCnt = BitUtil.pop(diff);
                int hetCnt = BitUtil.pop(hets);
                int sites = sameCnt + diffCnt - hetCnt;
                allDist[donor1][2][i]=(byte)diffCnt;
                allDist[donor1][3][i]=(byte)hetCnt;
                allDist[donor1][0][i]=(byte)sites;
                allDist[donor1][1][i]=(byte)sameCnt;
            }
        }
    }
    
    /**
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestInbredDonors(int targetTaxon, int startBlock, int endBlock, 
            int focusBlock, int[] donor1indices) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        int donorTaxaCnt=donorAlign.getSequenceCount();
        for (int d1 : donor1indices) {
            int testSites=0;
            int sameCnt = 0, diffCnt = 0, hetCnt = 0;
            for (int i = startBlock; i <=endBlock; i++) {
                sameCnt+=allDist[d1][1][i];
                diffCnt+=allDist[d1][2][i];
                hetCnt+=allDist[d1][3][i];
            }
            testSites= sameCnt + diffCnt - hetCnt;
            if(testSites<minTestSites) continue;     
            double testPropUnmatched = 1.0-(((double) (sameCnt) - (double)(0.5*hetCnt)) / (double) (testSites));
            int totalMendelianErrors=(int)((double)testSites*testPropUnmatched);
            inc+=1e-9;  //this looks strange, but makes the keys unique and ordered
            testPropUnmatched+=inc;  //this looks strange, but makes the keys unique and ordered 
            if(testPropUnmatched<lastKeytestPropUnmatched) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                    focusBlock, endBlock, testSites, totalMendelianErrors);
                DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                if(prev!=null) System.out.println("Inbred TreeMap index crash:"+testPropUnmatched);
                if(bestDonors.size()>maxDonorHypotheses) {
                    bestDonors.remove(bestDonors.lastKey());
                    lastKeytestPropUnmatched=bestDonors.lastKey();
                }
            }  
        }
        //TODO consider reranking by Poisson.  The calculating Poisson on the fly is too slow.
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
    }
    
    private int[] getAllBestDonorsAcrossChromosome(DonorHypoth[][] allDH, int minHypotheses) {
        TreeMap<Integer,Integer> bd=new TreeMap<Integer,Integer>();
        for (int i = 0; i < allDH.length; i++) {
            DonorHypoth dh=allDH[i][0];
            if(dh==null) continue;
            if(bd.containsKey(dh.donor1Taxon)) {
                bd.put(dh.donor1Taxon, bd.get(dh.donor1Taxon)+1);
            } else {
                bd.put(dh.donor1Taxon, 1);
            }
        }
        if(testing==1) System.out.println(bd.size()+":"+bd.toString());
        if(testing==1) System.out.println("");
        int highDonors=0;
        for (int i : bd.values()) {if(i>minHypotheses) highDonors++;}
        int[] result=new int[highDonors];    
        highDonors=0;
        for (Map.Entry<Integer,Integer> e: bd.entrySet()) {
            if(e.getValue()>minHypotheses) result[highDonors++]=e.getKey();}
        if(testing==1) System.out.println(Arrays.toString(result));
        return result;
    }
  
    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT, 
            int startBlock, int endBlock, int focusBlock, DonorHypoth[] inbredHypoth, boolean compareDonorDist) {
        int nonNull=0;
        for (DonorHypoth dH : inbredHypoth) {
            if(dH!=null) nonNull++;
        }
        int[] inbredIndices=new int[inbredHypoth.length];
        for (int i = 0; i < nonNull; i++) {
            inbredIndices[i]=inbredHypoth[i].donor1Taxon;   
        }
        return getBestHybridDonors(targetTaxon, mjT, mnT, startBlock,  endBlock, focusBlock, 
                inbredIndices, donorIndices, compareDonorDist);
    }
    
    /**
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT, 
            int startBlock, int endBlock, int focusBlock, int[] donor1Indices, int[] donor2Indices,
            boolean viterbiSearch) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        double[] donorDist={1.0,100.0};
       // Random r=new Random(1);
        for (int d1: donor1Indices) {
            long[] mj1=donorAlign.getAllelePresenceForSitesBlock(d1, 0, startBlock, endBlock+1);
            long[] mn1=donorAlign.getAllelePresenceForSitesBlock(d1, 1, startBlock, endBlock+1);
            for (int d2 : donor2Indices) {
                if((!viterbiSearch)&&(d1==d2)) continue;
                long[] mj2=donorAlign.getAllelePresenceForSitesBlock(d2, 0, startBlock, endBlock+1);
                long[] mn2=donorAlign.getAllelePresenceForSitesBlock(d2, 1, startBlock, endBlock+1);
                if(viterbiSearch) {
                    donorDist=IBSDistanceMatrix.computeHetBitDistances(mj1, mn1, mj2, mn2, minTestSites);
                    if((d1!=d2)&&(donorDist[0]<this.maximumInbredError)) continue;
                }
               // System.out.printf("%d %d %g %n",d1,d2,donorDist[0]);
    //            if(donorDist[0]<this.maximumInbredError) continue;  
                int[] mendErr=mendelErrorComparison(mjT, mnT, mj1, mn1, mj2, mn2);        
                if(mendErr[1]<minTestSites) continue;
                double testPropUnmatched=(double)(mendErr[0])/(double)mendErr[1];
                inc+=1e-9;
                testPropUnmatched+=inc; 
                if(testPropUnmatched<lastKeytestPropUnmatched) {
                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
                        focusBlock, endBlock, mendErr[1], mendErr[0]);
                    DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                    if(prev!=null) {
                        System.out.println("Hybrid TreeMap index crash:"+testPropUnmatched);
                    }
                    if(bestDonors.size()>maxDonorHypotheses) {
                        bestDonors.remove(bestDonors.lastKey());
                        lastKeytestPropUnmatched=bestDonors.lastKey();
                    }
                }  
            }
        }
        if(bestDonors.size()<1) {
            System.out.println("Houston problem here");
        }
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
    }
    
    private int[] mendelErrorComparison(long[] mjT, long[] mnT, long[] mj1, long[] mn1, 
            long[] mj2, long[] mn2) {
        int mjUnmatched=0;
        int mnUnmatched=0;
        int testSites=0;
        for (int i = 0; i < mjT.length; i++) {
            long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
            mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
            mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
            testSites+=Long.bitCount(siteMask);
        }
        int totalMendelianErrors=mjUnmatched+mnUnmatched;
       // double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
        return (new int[] {totalMendelianErrors, testSites});       
    }
    
        /**
     * Takes a donor hypothesis and applies it to the output alignment 
     * @param theDH
     * @param mna 
     */
    private void setAlignmentWithDonors(DonorHypoth[] theDH, boolean setJustFocus, MutableNucleotideAlignment mna) {
        if(theDH[0].targetTaxon<0) return;
        boolean print=false;
        int startSite=(setJustFocus)?theDH[0].getFocusStartSite():theDH[0].startSite;
        int endSite=(setJustFocus)?theDH[0].getFocusEndSite():theDH[0].endSite;
        if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
        if (print) System.out.println("B:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE; 
            for (int i = 0; (i < theDH.length) && (donorEst==Alignment.UNKNOWN_DIPLOID_ALLELE); i++) {
                if((theDH[i]==null)||(theDH[i].donor1Taxon<0)) continue;
                if(theDH[i].getErrorRate()>this.maximumInbredError) continue;
                byte bD1=donorAlign.getBase(theDH[i].donor1Taxon, cs);
                if(theDH[i].getPhaseForSite(cs)==0) {
                    donorEst=bD1;}
                else {
                    byte bD2=donorAlign.getBase(theDH[i].donor2Taxon, cs);
                    if(theDH[i].getPhaseForSite(cs)==2) {
                        donorEst=bD2;}
                    else if((bD1!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(bD2!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                        donorEst=(byte)((bD1&highMask)|(bD2&lowMask));
                    }
                }
            }
            byte knownBase=mna.getBase(theDH[0].targetTaxon, cs);
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                if(AlignmentUtils.isHeterozygous(donorEst)||AlignmentUtils.isHeterozygous(knownBase)) {
                    totalHets++;
                } else if(knownBase==donorEst) {
                    totalRight++;
                    siteCallCnt[cs]++;
                    taxonCallCnt[theDH[0].targetTaxon]++;
                } else {
                    totalWrong++;
  //                  System.out.printf("Known:%d donorEst: %d donorEst0:%d het:%s %n",knownBase,donorEst,donorEst0, AlignmentUtils.isHeterozygous(donorEst));
                    siteErrors[cs]++;
                    taxonErrors[theDH[0].targetTaxon]++;
                }
            }
            //TODO consider fixing the obvious errors.
            if(knownBase==Alignment.UNKNOWN_DIPLOID_ALLELE) mna.setBase(theDH[0].targetTaxon, cs, donorEst);
        } //end of cs loop
        if (print) System.out.println("E:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
    }
    
    private static void compareAlignment(String origFile, String maskFile, String impFile) {
        boolean taxaOut=false;
        Alignment oA=ImportUtils.readFromHapmap(origFile, false, (ProgressListener)null);
        System.out.printf("Orig taxa:%d sites:%d %n",oA.getSequenceCount(),oA.getSiteCount());        
        Alignment mA=ImportUtils.readFromHapmap(maskFile, false, (ProgressListener)null);
        System.out.printf("Mask taxa:%d sites:%d %n",mA.getSequenceCount(),mA.getSiteCount());
        Alignment iA=ImportUtils.readFromHapmap(impFile, false, (ProgressListener)null);
        System.out.printf("Imp taxa:%d sites:%d %n",iA.getSequenceCount(),iA.getSiteCount());
        int correct=0;
        int errors=0;
        int unimp=0;
        int hets=0;
        int gaps=0;
        for (int t = 0; t < iA.getSequenceCount(); t++) {
            int e=0,c=0,u=0,h=0;
            int oATaxa=oA.getIdGroup().whichIdNumber(iA.getFullTaxaName(t));
            for (int s = 0; s < iA.getSiteCount(); s++) {
                if(oA.getBase(oATaxa, s)!=mA.getBase(t, s)) {
                    byte ib=iA.getBase(t, s);
                    byte ob=oA.getBase(oATaxa, s);
                    if(ib==Alignment.UNKNOWN_DIPLOID_ALLELE) {unimp++; u++;}
                    else if(ib==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {gaps++;}
                    else if(ib==ob) {
                        correct++;
                        c++;
                    } else {
                        if(AlignmentUtils.isHeterozygous(ob)||AlignmentUtils.isHeterozygous(ib)) {hets++; h++;}
                        else {errors++; 
                            e++;
                            System.out.printf("%d %d %s %s %n",t,s,oA.getBaseAsString(oATaxa, s), iA.getBaseAsString(t, s));
                        }
                    }
                }       
            }
            if(taxaOut) System.out.printf("%s %d %d %d %d %n",iA.getTaxaName(t),u,h,c,e);
        }
        System.out.println("MFile\tIFile\tGap\tUnimp\tHets\tCorrect\tErrors");   
        System.out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%d%n",maskFile, impFile, gaps, unimp,hets,correct,errors);        
    }
      
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
      String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/impResults/";
      String rootIn="/Users/edbuckler/SolexaAnal/GBS/build20120701/impOrig/";
      
//      String root="/Volumes/LaCie/build20120701/impResults/";
//      String rootIn="/Volumes/LaCie/build20120701/impOrig/";

      
//        String origFile=rootIn+"all25.c10_s4096_s8191.hmp.txt.gz";
// //       String donorFile=rootIn+"w4096GapOf4_8KMerge20130507.hmp.txt.gz";
//        String donorFile=rootIn+"w4096Of4_8KMerge20130507.hmp.txt.gz";
//        String unImpTargetFile=rootIn+"all25.c10_s4096_s8191_masked.hmp.txt.gz";
//        String impTargetFile=root+"allInbredtest.c10_s4096_s8191.imp.hmp.txt.gz";
//        String impTargetFile2=root+"allHybridTest.c10_s4096_s8191.imp.hmp.txt.gz";
        

        String origFile=rootIn+"Z0CNtest.c10_s0_s24575.hmp.txt.gz";
        //String donorFile=rootIn+"w24575Of24KMerge20130513b.hmp.txt.gz";
        String donorFile=rootIn+"IMPw24575Of24KMerge20130513b.hmp.txt.gz";
        String unImpTargetFile=rootIn+"Z0CNtest.c10_s0_s24575_masked.hmp.txt.gz";
        String impTargetFile2=root+"HybridZ0CNtest.imp.hmp.txt.gz";
        
//        String origFile=rootIn+"10psample25.c10_s0_s24575.hmp.txt.gz";
//        String donorFile=rootIn+"IMPw24575Of24KMerge20130513b.hmp.txt.gz";
//        String unImpTargetFile=rootIn+"10psample25.c10_s0_s24575_masked.hmp.txt.gz";
//        String impTargetFile2=root+"10psample25.c10_s0_s24575_masked.imp.hmp.txt.gz";
        
//        String origFile=rootIn+"10psample25.c10_s0_s24575.hmp.txt.gz";
//        String donorFile=rootIn+"JustHMw24575OfHM224KMerge20130517b.hmp.txt.gz";
//        String unImpTargetFile=donorFile;
//        String impTargetFile2=root+"XJustHMw24575OfHM224KMerge20130517b.imp.hmp.txt.gz";
        
//        String origFile=rootIn+"SEEDsamp25.c10_s0_s24575.hmp.txt.gz";
//        String donorFile=rootIn+"IMPw24575Of24KMerge20130513b.hmp.txt.gz";
//        String unImpTargetFile=rootIn+"SEEDsamp25.c10_s0_s24575_masked.hmp.txt.gz";
//        String impTargetFile=root+"InbredSEED25.imp.hmp.txt.gz";
//        String impTargetFile2=root+"HybridSEED25.imp.hmp.txt.gz";
        
       
//      String origFile=rootIn+"10psample25.c10_s0_s24575.hmp.txt.gz";
//      String donorFile=rootIn+"Xw24575Of24KMerge20130513b.hmp.txt.gz";
//      String impTargetFile2=root+"IMPw24575Of24KMerge20130513b.hmp.txt.gz";

//        
//        String origFile=rootIn+"ABQTL25.c10_s0_s24575.hmp.txt.gz";
//        String donorFile=rootIn+"w24575Of24KMerge20130513b.hmp.txt.gz";
//        String unImpTargetFile=rootIn+"ABQTL25.c10_s0_s24575_masked.hmp.txt.gz";
//        String impTargetFile=root+"ABQTLtest.imp.hmp.txt.gz";
//        String impTargetFile2=root+"HybridABQTLtest.imp.hmp.txt.gz";
      
//        String origFile=rootIn+"Z0CN26.c10_s0_s24575.hmp.txt.gz";
////        String donorFile=rootIn+"JustHMw24575OfHM224KMerge20130517b.imp.hmp.txt.gz";
//        String donorFile=rootIn+"JustHMw24575OfHM224KMerge20130517b.hmp.txt.gz";
//        String unImpTargetFile=rootIn+"Z0CN26.c10_s0_s24575.hmp.txt.gz";
//        String impTargetFile2=root+"Z0CN26.c10_s0_s24575.imp.hmp.txt.gz";
        
//        String origFile=rootIn+"10psample26.c10_s0_s24575.hmp.txt.gz";
//        String donorFile=rootIn+"JustHMw24575OfHM224KMerge20130517b.imp.hmp.txt.gz";
////        String donorFile=rootIn+"JustHMw24575OfHM224KMerge20130517b.hmp.txt.gz";
//        String unImpTargetFile=rootIn+"10psample26.c10_s0_s24575.hmp.txt.gz";
//        String impTargetFile2=root+"10psample26.c10_s0_s24575.imp.hmp.txt.gz";


        MinorWindowViterbiImputation e64NNI=new MinorWindowViterbiImputation(donorFile, unImpTargetFile, impTargetFile2, 20, 50, 100, 0.01, true, false);
        compareAlignment(origFile,unImpTargetFile,impTargetFile2);
        
//        String origFile=rootIn+"10psample25.c10_s0_s24575.hmp.txt.gz";
//        String donorFile=rootIn+"JustHMw24575OfHM224KMerge20130517b.hmp.txt.gz";
//        String unImpTargetFile=donorFile;
//        String impTargetFile2=root+"JustHMw24575OfHM224KMerge20130517b.imp.hmp.txt.gz";
//
//        MinorWindowViterbiImputation e64NNI=new MinorWindowViterbiImputation(donorFile, unImpTargetFile, impTargetFile2, 20, 50, 100, 0.01, true, true);
//        
        System.out.println("Resolve Method 0: Minor 20");
        
//        e64NNI=new MinorWindowViterbiImputation(donorFile, unImpTargetFile, impTargetFile, 20, false);
//        compareAlignment(origFile,unImpTargetFile,impTargetFile);
        
//        e64NNI=new MinorWindowViterbiImputation(donorFile, impTargetFile, impTargetFile2, 20, true);
//        compareAlignment(origFile,unImpTargetFile,impTargetFile2);
        for (double d = 0.01; d < 0.02; d+=0.01) {
            System.out.println("MaxHybrid Error Rate:"+d);
            
            
        }
//        e64NNI=new MinorWindowViterbiImputation(donorFile, unImpTargetFile, impTargetFile2, 20, 100, true);
//        compareAlignment(origFile,unImpTargetFile,impTargetFile2);
        


    }
    
}
