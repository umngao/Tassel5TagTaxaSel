/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.dna.tag.Tag;

/**
 *  Class needed for storing reference tags into RepGen SQLite
 *  tables
 * @author lcj34
 *
 */
public class RefTagData {
    private final Tag myTag;
    private final String chromosome;
    private final int position;
    private final int refGenomeID;

    public RefTagData(Tag myTag, String chromosome,int position, int refGenomeID) {
        this.myTag = myTag;
        this.chromosome = chromosome;
        this.position = position;
        this.refGenomeID = refGenomeID;
    }

    public Tag tag() {
        return myTag;
    }

    public String chromosome() {
        return chromosome;
    }
    public int position() {
        return position;
    }
    public int refGenomeID() {
        return refGenomeID;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        RefTagData that = (RefTagData) obj;

        if (!(tag().equals(that.tag()))) return false;
        if (!(chromosome().equals(that.chromosome()))) return false;
        if (position() != that.position()) return false;
        if (refGenomeID() != that.refGenomeID()) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 37 * hash + this.myTag.hashCode();
        try {
            hash = 37 * hash + Integer.parseInt(chromosome);
        } catch (NumberFormatException nfe) {
            // add nothing for chromsome if not int
        }
        hash = 37 * hash + this.position;
        hash = 37 * hash + this.refGenomeID;
        return hash;
    }

}
