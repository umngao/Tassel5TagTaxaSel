/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.read;

import net.maizegenetics.dna.BaseEncoder;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 * Holding paired end read. Providing merging function
 * @author Fei Lu
 */
public class PERead {
    Read rf;
    Read rb;
    Read pContig = null;
    
    public PERead (Read rf, Read rb) {
        this.rf = rf;
        this.rb = rb;
    }
    
    public boolean merge (boolean ifPrintAlignment) {
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 10, (short) 10);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        int minOverlap = 10;
        double minIden = 0.6;
        String queryS = this.rf.getSeq();
        DNASequence query = new DNASequence(queryS);
        String hitS = this.rb.getReverseComplementarySeq();
        String hitQualS = this.rb.getReverseQual();
        DNASequence hit = new DNASequence(hitS);
        //int halfLength = queryS.length()/2;
        SequencePair<DNASequence, NucleotideCompound> psa = null;
        psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
        int queryStart=0;
        try {
            queryStart = psa.getIndexInQueryAt(1);
        }
        catch (NullPointerException e) {
            //When PE has bad sequences, e.g. NNNNNNNNNNNNNN. No alignment in pas
            return false;
        }
        int queryEnd = psa.getIndexInQueryAt(psa.getLength());
        int hitStart = psa.getIndexInTargetAt(1);
        int hitEnd = psa.getIndexInTargetAt(psa.getLength());
        int overlap = psa.getLength();
        int idenNum = psa.getNumIdenticals();
        if (hitStart > 5) return false;
        if (queryEnd < queryS.length()-5) return false;
        if (overlap < minOverlap) return false;
        double iden = (double)idenNum/overlap;
        //if (iden < minIden) return false;
        StringBuilder sbSeq = new StringBuilder();
        StringBuilder sbQual = new StringBuilder();
        sbSeq.append(queryS.substring(0, queryEnd));
        sbQual.append(this.rf.getQual().substring(0, queryEnd));
        sbSeq.append(hitS.substring(hitEnd-1, hitS.length()));
        sbQual.append(hitQualS.substring(hitEnd-1, hitS.length()));
        String contigS = sbSeq.toString();
        String contigQualS = sbQual.toString();
        String ID = "@"+this.rf.ID.replaceFirst("@", "") + this.rb.ID.replaceFirst("@", "_____");
        String des = "+"+this.rf.des.replaceFirst("\\+", "") + this.rb.des.replaceFirst("\\+", "_____");
        pContig = new Read(ID, contigS, des, contigQualS);
        if (ifPrintAlignment) {
            System.out.println("********************************************************\n");
            System.out.println("QueryStart:\t" + String.valueOf(queryStart));
            System.out.println("QueryEnd:\t" + String.valueOf(queryEnd));
            System.out.println("QueryLength:\t" + String.valueOf(rf.getReadLength()));
            System.out.println("HitStart:\t" + String.valueOf(hitStart));
            System.out.println("HitEnd:\t" + String.valueOf(hitEnd));
            System.out.println("HitLength:\t" + String.valueOf(rb.getReadLength()));
            System.out.println("OverlapLength:\t" + String.valueOf(psa.getLength()));
            System.out.println("Identity:\t" + String.valueOf(iden)+"\t" + String.valueOf(idenNum));
            System.out.println("PE contig length:\t" + String.valueOf(contigS.length()));
            System.out.println(psa.toString(1000));
            System.out.println(contigS);
            System.out.println(queryS);
            System.out.println(hitS+"\n");
            System.out.println(contigQualS);
            System.out.println(this.rf.getQual());
            System.out.println(hitQualS);
            System.out.println("\n\n");
        }
        return true;
    }
    
    /**
     * Return forward read
     * @return 
     */
    public Read getForwardRead () {
        return rf;
    }
    
    /**
     * Return backward read
     * @return 
     */
    public Read getBackwardRead () {
        return rb;
    }
    
    /**
     * Return Pcontig read
     * Return null when there is no Pcontig
     * @return 
     */
    public Read getPContig () {
        return pContig;
    }
}
