package net.maizegenetics.dna.tag;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.dna.BaseEncoder;

import java.nio.ByteBuffer;

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
}
