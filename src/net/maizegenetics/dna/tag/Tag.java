package net.maizegenetics.dna.tag;

import net.maizegenetics.dna.BaseEncoder;

import java.io.Serializable;

/**
 * Interface for tags that captures these bit encoded sequence and there length.
 * 
 * @author Ed Buckler
 */
public interface Tag {

    public String sequence();

    public long[] seq2Bit();

    public byte[] seq2BitAsBytes();

    public short seqLength();

    public boolean isReference();

    default String toCSVString() {
        return sequence() + "," + seqLength();
    }

    default String toReverseComplement() {
        return BaseEncoder.getReverseComplement(sequence());
    }

}
