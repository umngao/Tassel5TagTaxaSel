/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.dna.tag.Tag;

/**
 * @author lcj34
 *
 */
public class TagCorrelationInfo {

    private final Tag tag2;
    private final double t1t2_pearson;
    private final double t1t2_spearman;
    private final double t1pt2p_pearson;
    private final double t1pt2p_r2;
    
    public TagCorrelationInfo(Tag tag2, double t1t2_pearson, double t1t2_spearman, double t1pt2p_pearson, double r2) {
        this.tag2 = tag2;
        this.t1t2_pearson = t1t2_pearson;
        this.t1t2_spearman = t1t2_spearman;
        this.t1pt2p_pearson = t1pt2p_pearson;
        this.t1pt2p_r2 = r2;
    }
    
    public Tag getTag2() {
        return tag2;
    }
    
    public double getT1t2_pearson() {
        return t1t2_pearson;
    }
    public double getT1t2_spearman() {
        return t1t2_spearman;
    }
    public double getT1pt2p_pearson() {
        return t1pt2p_pearson;
    }
    public double r2() {
        return t1pt2p_r2;
    }
}
