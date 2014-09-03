package net.maizegenetics.dna.tag;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.dna.BaseEncoder;

import java.io.Serializable;
import java.nio.ByteBuffer;
import java.util.Arrays;

/**
 * Builder for tags that optimizes the memory footprint.
 * Tags are encoded in long with 2 bits per bp.  In GBS, we generally use only record the 64 or 96bp.
 *
 * Custom classes are generated with 32, 64 or 96bp tags (1, 2 or 3 longs) otherwise an array is used.
 *
 * @author Ed Buckler
 */
public class TagBuilder {
    private TagBuilder() {}

    public static Tag instance(long[] seq2Bit, short length) {
        switch (seq2Bit.length) {
            case 0: return null;
            case 1: return new Tag1Long(seq2Bit,(byte)length);
            case 2: return new Tag2Long(seq2Bit,(byte)length);
            case 3: return new Tag3Long(seq2Bit,(byte)length);
            default: return new TagVarLong(seq2Bit,length);
        }
    }

    public static Tag instance(byte[] seq2BitInBytes, short length) {
        int seqBitLength=seq2BitInBytes.length/8;
        long[] seq2Bit=new long[seqBitLength];
        ByteBuffer bb=ByteBuffer.wrap(seq2BitInBytes);
        for (int i = 0; i < seq2Bit.length; i++) {
            seq2Bit[i]=bb.getLong();
        }
        return instance(seq2Bit,length);
    }

    public static Tag instance(String sequence) {
        long[] seq2Bit= AbstractTag.getLongArrayFromSeq(sequence);
        if(seq2Bit==null) return null;
        int length=sequence.length();
        return instance(seq2Bit,(short)length);
    }
}

class Tag1Long extends AbstractTag {
    //memory 8 + 8 + 1 =17 bytes
    //An array would add 12 to it
    private final long val0;
    private final byte length;

    Tag1Long(long[] val, byte length) {
        val0=val[0];
        this.length=length;
    }

    @Override
    public long[] seq2Bit() {
        return new long[]{val0};
    }

    @Override
    public short seqLength() {
        return length;
    }
}

class Tag2Long extends AbstractTag {
    //memory 8 + 16 + 1 =25 bytes
    //An array would add 12 to it
    private final long val0, val1;
    private final byte length;

    Tag2Long(long[] val, byte length) {
        val0=val[0];
        val1=val[1];
        this.length=length;
    }

    @Override
    public long[] seq2Bit() {
        return new long[]{val0,val1};
    }

    @Override
    public short seqLength() {
        return length;
    }
}

class Tag3Long extends AbstractTag {
    //memory 8 + 24 + 1 = 33 bytes
    //An array would add 12 to it
    private long val0, val1, val2;
    private byte length;

    Tag3Long(long[] val, byte length) {
        val0=val[0];
        val1=val[1];
        val2=val[2];
        this.length=length;
    }

    @Override
    public long[] seq2Bit() {
        return new long[]{val0,val1,val2};
    }

    @Override
    public short seqLength() {
        return length;
    }

}

class TagVarLong extends AbstractTag {
    //memory 8 + 12 + 8*LongLen + 2 = XX bytes
    private long[] val;
    private short length;

    TagVarLong(long[] val, short length) {
        this.val=val;
        this.length=length;
    }

    @Override
    public long[] seq2Bit() {
        return val;
    }

    @Override
    public short seqLength() {
        return length;
    }

}
