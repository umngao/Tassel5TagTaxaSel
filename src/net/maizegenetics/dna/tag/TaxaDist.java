package net.maizegenetics.dna.tag;

import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Multiset;
import com.google.common.primitives.Shorts;


/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
public class TaxaDist {
    //minimal size 8 + 12 + 12 + 10 + 12 + 4+ 4 + 4 = 66
    //This could be changed for the singletons by just making a new class
    private short[][] taxaWithTag;
    private int[] size;
    private int totalSize;
    private final int maxTaxa;
    private int singleTaxa=Integer.MIN_VALUE;  //this halves the memory footprint for the singleton tags

    public TaxaDist(int maxTaxa) {
        this.maxTaxa=maxTaxa;
        initializeArrays();
    }

    public TaxaDist(int maxTaxa, int taxaWithTag) {
        this.maxTaxa=maxTaxa;
        singleTaxa=maxTaxa;
        totalSize=1;
    }

    public TaxaDist(TaxaDist td) {
        this(td.maxTaxa);
        for (Integer tn: td.taxaDepthMap()) {
            increment(tn);
        }
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

    public synchronized TaxaDist increment(int taxaNum) {
        int shortSet=taxaNum>>>16;
        taxaWithTag[shortSet]=Shorts.ensureCapacity(taxaWithTag[shortSet],size[shortSet]+1,size[shortSet]);
        taxaWithTag[shortSet][size[shortSet]]=(short)taxaNum;
        size[shortSet]++;
        totalSize++;
        return this;
    }

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

    public int[] encodeTaxaDepth() {
        int[][] tds=taxaWithDepths();
        int[] result=new int[tds[0].length];
        for (int i = 0; i < tds[0].length; i++) {
            result[i]=(tds[0][i]<<8)|(tds[1][i]);
        }
        return result;
    }

    public Multiset<Integer> taxaDepthMap() {
        int[][] tds=taxaWithDepths();
        ImmutableMultiset.Builder<Integer> taxaCnts = new ImmutableMultiset.Builder<>();
        for (int i = 0; i < tds[0].length; i++) {
            taxaCnts.setCount(tds[0][i],tds[1][i]);
        }
        return taxaCnts.build();
    }

    public int totalDepth(){
        return totalSize;
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
