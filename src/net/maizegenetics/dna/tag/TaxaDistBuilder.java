package net.maizegenetics.dna.tag;

import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Multiset;
import com.google.common.primitives.Shorts;


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

}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistExpandable implements TaxaDistribution {
    //minimal size 8 + 12 + 12 + 10 + 12 + 4+ 4 + 4 = 66
    //This could be changed for the singletons by just making a new class
    private short[][] taxaWithTag;
    private int[] size;
    private int totalSize;
    private final int maxTaxa;
    private int singleTaxa=Integer.MIN_VALUE;  //this halves the memory footprint for the singleton tags

    public TaxaDistExpandable(int maxTaxa) {
        this.maxTaxa=maxTaxa;
        initializeArrays();
    }

    public TaxaDistExpandable(int maxTaxa, int taxaWithTag) {
        this.maxTaxa=maxTaxa;
        singleTaxa=maxTaxa;
        totalSize=1;
    }

    private void initializeArrays() {
        int shortSets=1+(maxTaxa>>>16);
        taxaWithTag=new short[shortSets][];
        for (int i = 0; i < shortSets; i++) {
            taxaWithTag[i]=new short[5];
        }
        size=new int[shortSets];
        if(singleTaxa>=0) increment(singleTaxa);
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
        if(totalSize==1) {depths[singleTaxa]=1; return depths;}
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
        int[] result=new int[tds[0].length];
        for (int i = 0; i < tds[0].length; i++) {
            result[i]=(tds[0][i]<<8)|(tds[1][i]);
        }
        return result;
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
    public String toString() {
        return "TaxaDist{" +
                "totalSize=" + 1 +
                ", maxTaxa=" + maxTaxa +
                ", "+taxaDepthMap().toString()+
                '}';
    }
}