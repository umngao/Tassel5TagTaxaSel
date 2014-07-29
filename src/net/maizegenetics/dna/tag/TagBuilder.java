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
 * Custom classes are generated with 64 or 96bp tags (2 or 3 longs) otherwise an array is used.
 *
 * @author Ed Buckler
 */
public class TagBuilder {
    private TagBuilder() {}

    public static Tag instance(long[] seq2Bit, short length) {
        if(seq2Bit.length==2) return new Tag2Long(seq2Bit,(byte)length);
        if(seq2Bit.length==3) return new Tag3Long(seq2Bit,(byte)length);
        if(seq2Bit.length>3) return new TagVarLong(seq2Bit,(short)length);
        return null;
    }

    public static Tag instance(byte[] seq2BitInBytes, short length) {
        int seqBitLength=seq2BitInBytes.length/8;
        long[] seq2Bit=new long[seqBitLength];
        ByteBuffer bb=ByteBuffer.wrap(seq2BitInBytes);
        for (int i = 0; i < seq2Bit.length; i++) {
            seq2Bit[i]=bb.getLong();
        }
        if(seq2Bit.length==2) return new Tag2Long(seq2Bit,(byte)length);
        if(seq2Bit.length==3) return new Tag3Long(seq2Bit,(byte)length);
        if(seq2Bit.length>3) return new TagVarLong(seq2Bit,(short)length);
        return null;
    }

    public static Tag instance(String sequence) {
        long[] seq2Bit= BaseEncoder.getLongArrayFromSeq(sequence);
        int length=sequence.length();
        if(seq2Bit.length==2) return new Tag2Long(seq2Bit,(byte)length);
        if(seq2Bit.length==3) return new Tag3Long(seq2Bit,(byte)length);
        return new TagVarLong(seq2Bit,(short)length);
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Tag2Long tag2Long = (Tag2Long) o;

        if (length != tag2Long.length) return false;
        if (val0 != tag2Long.val0) return false;
        if (val1 != tag2Long.val1) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = (int) (val0 ^ (val0 >>> 32));
        result = 31 * result + (int) (val1 ^ (val1 >>> 32));
        result = 31 * result + (int) length;
        return result;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Tag3Long tag3Long = (Tag3Long) o;

        if (length != tag3Long.length) return false;
        if (val0 != tag3Long.val0) return false;
        if (val1 != tag3Long.val1) return false;
        if (val2 != tag3Long.val2) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = (int) (val0 ^ (val0 >>> 32));
        result = 31 * result + (int) (val1 ^ (val1 >>> 32));
        result = 31 * result + (int) (val2 ^ (val2 >>> 32));
        result = 31 * result + (int) length;
        return result;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        TagVarLong that = (TagVarLong) o;

        if (length != that.length) return false;
        if (!Arrays.equals(val, that.val)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(val);
        result = 31 * result + (int) length;
        return result;
    }
}
