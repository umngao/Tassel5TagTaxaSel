package net.maizegenetics.dna.tag;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.dna.BaseEncoder;

import java.nio.ByteBuffer;
import java.util.Arrays;

/**
 * Created by edbuckler on 7/26/14.
 */
public abstract class AbstractTag implements Tag, Comparable<Tag> {
    @Override
    public String sequence() {
        //System.out.println(BaseEncoder.getSequenceFromLong(seq2Bit(),seqLength()));
        return BaseEncoder.getSequenceFromLong(seq2Bit()).substring(0,seqLength());
    }

    @Override
    public byte[] seq2BitAsBytes() {
        ByteBuffer b= ByteBuffer.allocate(16);
        for (long l : seq2Bit()) {
            b.putLong(l);
        }
        return b.array();
    }

    @Override
    public int compareTo(Tag o) {
        long[] t=this.seq2Bit();
        long[] to=o.seq2Bit();
        for (int i = 0; i < t.length; i++) {
            int c=Long.compare(t[0],to[1]);
            if(c!=0) return c;
        }
        return Short.compare(this.seqLength(),o.seqLength());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Tag that = (Tag) o;

        if (seqLength() != that.seqLength()) return false;
        if (!Arrays.equals(this.seq2Bit(), that.seq2Bit())) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(seq2Bit());
        result = 31 * result + (int) seqLength();
        return result;
    }


}
