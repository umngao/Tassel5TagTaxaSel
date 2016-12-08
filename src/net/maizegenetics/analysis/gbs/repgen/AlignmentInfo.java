/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.dna.tag.Tag;

/**
 * This class used by RepGenAlignerPlugin to store alignment
 * info to db table tagAlignments
 * 
 * @author lcj34
 *
 */
public class AlignmentInfo implements Comparable<AlignmentInfo>{
    private final Tag tag2;
    private  final String tag2chrom;
    private  final int tag2pos;
    private  final int alignmentPos;
    private final int ref_strand; // forward/plus=1, reverse/minus=0, unknown = -1,(mostly) as per Position interface
    private  final int myScore;

    public AlignmentInfo(Tag tag2, String chromosome, int position, int alignmentpos, int ref_strand, int score) {
        this.tag2 = tag2;
        this.tag2chrom = chromosome;
        this.tag2pos = position;
        this.alignmentPos = alignmentpos;
        this.ref_strand = ref_strand;
        this.myScore = score;
    }

    public Tag tag2() {
        return tag2;
    }
    public String tag2chrom() {
        return tag2chrom;
    }

    public int tag2pos() {
        return tag2pos;
    }
    
    public int alignmentPos() {
        return alignmentPos;
    }
    public int ref_strand() {
        return ref_strand;
    }
    public  int score() {
        return myScore;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Alignment:");
        sb.append("\tTag2:").append(tag2.sequence());
        sb.append("\tChr:").append(tag2chrom);
        sb.append("\tPos:").append(tag2pos);
        sb.append("\tAlignmentPos:").append(alignmentPos);
        sb.append("\tScore:").append(myScore);
        sb.append("\n");
        return sb.toString();
    }
    @Override
    public int compareTo(AlignmentInfo other) {
        // TODO Auto-generated method stub
        return this.myScore > other.score() ? 1: this.myScore < other.score() ? -1 : 0;
    }
}
