package net.maizegenetics.dna.tag;

import java.io.Serializable;

/**
 * Interface for tags that captures these bit encoded sequence and there length.
 * 
 * @author Ed Buckler
 */
public interface Tag {

    public long[] seq2Bit();

    public short seqLength();

}
