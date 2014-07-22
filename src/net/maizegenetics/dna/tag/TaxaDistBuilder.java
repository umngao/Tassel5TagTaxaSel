package net.maizegenetics.dna.tag;

import cern.colt.list.IntArrayList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Multiset;
import com.google.common.primitives.Shorts;
import com.google.common.primitives.UnsignedBytes;

import java.io.Serializable;


/**
 * Builder for TaxaDistribution.  Deals with the
 */
public class TaxaDistBuilder {
    private TaxaDistBuilder() {}

    /**
     * Create a TaxaDistribution with only a single taxa with a tag.  Very memory efficient, but needs conversion add
     * an additional taxon with a tag
     * @param maxTaxa
     * @param taxaWithTag
     * @return
     */
    public static TaxaDistribution create(int maxTaxa, int taxaWithTag) {
        return new TaxaDistSingleTaxon(maxTaxa,taxaWithTag);
    }

    /**
     * Create an expandable TaxaDistribution with no initial values
     * @param maxTaxa
     * @return
     */
    public static TaxaDistribution create(int maxTaxa) {
        return new TaxaDistExpandable(maxTaxa);

    }

    /**
     * Copies a TaxaDistribution to an expandable TaxaDistribution.  Can be used to convert a single TaxaDistribution
     * to an expandable version
     * @param srcTaxaDist
     * @return
     */
    public static TaxaDistribution create(TaxaDistribution srcTaxaDist) {
        TaxaDistribution dstTD=new TaxaDistExpandable(srcTaxaDist.maxTaxa());
        int[] depths=srcTaxaDist.depths();
        for (int taxaIndex = 0; taxaIndex < depths.length; taxaIndex++) {
            for (int j = 0; j < depths[taxaIndex]; j++) {
                dstTD.increment(taxaIndex);
            }
        }
        return dstTD;
    }

    /**
     * Create expandable TaxaDistribution from encoded taxa distribution.  The int[] encoding use first three bytes
     * for taxa index, and last byte for depth as unsigned byte.  Depths greater than 256 just increase.
     * @param maxTaxa
     * @param encodeTaxaDepths
     * @return
     */
    public static TaxaDistribution create(int maxTaxa, int[] encodeTaxaDepths) {
        TaxaDistribution dstTD=new TaxaDistExpandable(maxTaxa);
        for (int taxaDepth : encodeTaxaDepths) {
            int taxa=taxaDepth>>>8;
            int depth=UnsignedBytes.toInt((byte)taxaDepth);
            for (int i = 0; i < depth; i++) {
                dstTD.increment(taxa);
            }
        }
        return dstTD;
    }

}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistExpandable implements TaxaDistribution, Serializable {
    //minimal size 8 + 12 + 12 + 10 + 12 + 4+ 4 = 66
    //This could be changed for the singletons by just making a new class
    private short[][] taxaWithTag;
    private int[] size;
    private int totalSize;
    private final int maxTaxa;

    public TaxaDistExpandable(int maxTaxa) {
        this.maxTaxa=maxTaxa;
        initializeArrays();
    }

    private void initializeArrays() {
        int shortSets=1+(maxTaxa>>>16);
        taxaWithTag=new short[shortSets][];
        for (int i = 0; i < shortSets; i++) {
            taxaWithTag[i]=new short[5];
        }
        size=new int[shortSets];
    }

    @Override
    public synchronized TaxaDistribution increment(int taxaNum) {
        int shortSet=taxaNum>>>16;
        taxaWithTag[shortSet]=Shorts.ensureCapacity(taxaWithTag[shortSet],size[shortSet]+1,size[shortSet]);
        taxaWithTag[shortSet][size[shortSet]]=(short)taxaNum;
        size[shortSet]++;
        totalSize++;
        return this;
    }

    @Override
    public int[] depths() {
        int[] depths=new int[maxTaxa];
        for (int i = 0; i < taxaWithTag.length; i++) {
            int base=i*(1<<16);
            for (int j = 0; j < size[i]; j++) {
                depths[unSignShort(taxaWithTag[i][j])+base]++;
            }
        }
        return depths;
    }

    @Override
    public int[][] taxaWithDepths() {
        int[] depths=depths();
        int countNoZero=0;
        for (int depth : depths) if(depth>0) countNoZero++;
        int[][] taxaDepth=new int[2][countNoZero];
        countNoZero=0;
        for (int i = 0; i <depths.length ; i++) {
            if(depths[i]>0) {
                taxaDepth[0][countNoZero]=i;
                taxaDepth[1][countNoZero]=depths[i];
                countNoZero++;
            }
        }
        return taxaDepth;
    }

    @Override
    public int[] encodeTaxaDepth() {
        int[][] tds=taxaWithDepths();
        IntArrayList result=new IntArrayList(tds[0].length);
        for (int i = 0; i < tds[0].length; i++) {
            int depth=tds[1][i];
            while(depth>0) {
                byte outDepth=(depth<256)?UnsignedBytes.checkedCast(depth):UnsignedBytes.checkedCast(255);
                result.add((tds[0][i]<<8)|(outDepth));
                depth-=255;
            }
        }
        result.trimToSize();
        return result.elements();
    }

    @Override
    public Multiset<Integer> taxaDepthMap() {
        int[][] tds=taxaWithDepths();
        ImmutableMultiset.Builder<Integer> taxaCnts = new ImmutableMultiset.Builder<>();
        for (int i = 0; i < tds[0].length; i++) {
            taxaCnts.setCount(tds[0][i],tds[1][i]);
        }
        return taxaCnts.build();
    }

    @Override
    public int totalDepth(){
        return totalSize;
    }

    @Override
    public int maxTaxa() {
        return maxTaxa;
    }

    @Override
    public int memorySize() {
        //minimal size 8 (object) + 12 (outer short array) + 12 (sizeArray) + 4+ 4 = 40
        int size=40;
        for (int i = 0; i < taxaWithTag.length; i++) {
            size+=16+(taxaWithTag[i].length*2);  //4 size array + 12 inner short array plus of the size of it
        }
        return size;
    }

    private int unSignShort(short v) {
        if(v<0) return -v+1;
        return v;
    }

    @Override
    public String toString() {
        return "TaxaDist{" +
                "totalSize=" + totalSize +
                ", maxTaxa=" + maxTaxa +
                ", "+taxaDepthMap().toString()+
                '}';
    }
}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistSingleTaxon implements TaxaDistribution {
    //minimal size 8 + 4 + 4  = 16
    private final int maxTaxa;
    private final int singleTaxa;


    public TaxaDistSingleTaxon(int maxTaxa, int taxaWithTag) {
        this.maxTaxa=maxTaxa;
        singleTaxa=taxaWithTag;
    }

    @Override
    public synchronized TaxaDistribution increment(int taxaNum) {
        throw new UnsupportedOperationException("TaxaDistSingleTaxon cannot be increment.  Expand first.");
    }

    @Override
    public int[] depths() {
        int[] depths=new int[maxTaxa];
        depths[singleTaxa]=1;
        return depths;
    }

    @Override
    public int[][] taxaWithDepths() {
        return new int[][]{{singleTaxa},{1}};
    }

    @Override
    public int[] encodeTaxaDepth() {
        int[][] tds=taxaWithDepths();
        int[] result=new int[tds[0].length];
        for (int i = 0; i < tds[0].length; i++) {
            result[i]=(tds[0][i]<<8)|(tds[1][i]);
        }
        return result;
    }

    @Override
    public Multiset<Integer> taxaDepthMap() {
        return new ImmutableMultiset.Builder<Integer>().add(singleTaxa).build();
    }

    @Override
    public int totalDepth(){
        return 1;
    }

    @Override
    public int maxTaxa() {
        return maxTaxa;
    }

    @Override
    public int memorySize() {
        return 16;
    }

    @Override
    public String toString() {
        return "TaxaDist{" +
                "totalSize=" + 1 +
                ", maxTaxa=" + maxTaxa +
                ", "+taxaDepthMap().toString()+
                '}';
    }
}