/*
 * TagsByTaxaByteHDF5TaxaGroups
 */
package net.maizegenetics.dna.tag;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.BiMap;
import net.maizegenetics.util.*;

import java.util.concurrent.ExecutionException;

/**
 * Tags by Taxa file based on the HDF5 data structure.  It can work with HDF5 files with either chunk directions.
 * 
 * @author Ed Buckler
 */
public class TagsByTaxaHDF5 extends AbstractTagsByTaxa {
    public enum Chunking {
        Taxa, Tags
    }
    private final boolean inTaxaChunking;
    private final int tagCount;
    private final IHDF5Reader reader;
    private final BiMap<String, Integer> taxaNameIndexBiMap;


    //todo TAS-203 this needs to be redone for chunks for taxa, as at 98M tags this will never get pivoted.
    private LoadingCache<Long, byte[]> taxaCache=null;

    private CacheLoader<Long, byte[]> taxaLoader = new CacheLoader<Long, byte[]>() {
        public byte[] load(Long key) {
           return HDF5Utils.getTagDistBlockForTaxon(reader, getTaxonFromKey(key), getTaxaBlockFromKey(key)); //HDF5Utils
        }
    };

    private LoadingCache<Integer, byte[]> tagCache=null;

    private CacheLoader<Integer, byte[]> tagLoader = new CacheLoader<Integer, byte[]>() {
        public byte[] load(Integer key) {
            return HDF5Utils.getTaxaDistForTag(reader, key); //HDF5Utils
        }
    };

    public TagsByTaxaHDF5(String infile) {
        reader = HDF5Factory.openForReading(infile);
        tagCount = HDF5Utils.getHDF5TagCount(reader);
        inTaxaChunking= HDF5Utils.isTagsByTaxaInTaxaDirection(reader);
        tagLengthInLong = HDF5Utils.getHDF5TagLengthInLong(reader);
        taxaNum = HDF5Utils.getNumberOfTaxaInTBT(reader);
        this.tags = HDF5Utils.getTags(reader);
        taxaNameIndexBiMap=HDF5Utils.getTBTMapOfRowIndices(reader);
        taxaList=HDF5Utils.getTaxaListInTBTOrder(reader);
        //todo need to make tag lengths an integer
        int[] tL=HDF5Utils.getTagLengths(reader);
        this.tagLength=new byte[tL.length];
        for (int i = 0; i < tagLength.length; i++) {
            tagLength[i]=(byte)tL[i];
        }

        if(inTaxaChunking) {
            taxaCache = CacheBuilder.newBuilder()
                    .maximumSize(1000)
                    .build(taxaLoader);
        } else {
            tagCache = CacheBuilder.newBuilder()
                    .maximumSize(100_000)
                    .build(tagLoader);
        }
    }


    private long getTaxaCacheKey(int taxonIndex, int tagIndex) {
        return ((long) taxonIndex << 32) | (long) (tagIndex/Tassel5HDF5Constants.BLOCK_SIZE);
    }

    private int getTaxonFromKey(long key) {
        return (int) (key >>> 32);
    }

    private int getTaxaBlockFromKey(long key) {
        return (int) ((key << 32) >>> 32);
    }

    private int getOffsetWithinBlock(int tagIndex) {
        return tagIndex%Tassel5HDF5Constants.BLOCK_SIZE;
    }


    @Override
    public int getIndexOfTaxaName(String taxon) {
        return taxaNameIndexBiMap.get(taxon);
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        if(inTaxaChunking) {
            return getReadForTaxon(taxaIndex,tagIndex);
        } else {
            return getTaxaReadCountsForTag(tagIndex)[taxaIndex];
        }

    }

    public byte[] getReadCountDistributionForTaxon(int taxonIndex) {
        byte[] readsCnts=new byte[tagCount];
        for (int i = 0; i < readsCnts.length; i++) {
            readsCnts[i]=getReadForTaxon(i,taxonIndex);
        }
        return readsCnts;
    }


    private byte getReadForTaxon(int tagIndex, int taxonIndex) {
        try {
            return taxaCache.get(getTaxaCacheKey(taxonIndex,tagIndex))[getOffsetWithinBlock(tagIndex)];
        } catch (ExecutionException e) {
            e.printStackTrace();
            return -1;
        }
    }

    @Override
    public byte[] getTaxaReadCountsForTag(int tagIndex) {
        try {
            return tagCache.get(tagIndex);
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public void setMethodByRows(boolean rowSetMethod) {
        throw new UnsupportedOperationException("Not supported");
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
            throw new UnsupportedOperationException("Not supported");
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        throw new UnsupportedOperationException("Not supported");
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        throw new UnsupportedOperationException("Not supported");
    }

    @Override
    public void getFileReadyForClosing() {
            throw new UnsupportedOperationException("Not supported");
    }

    @Override
    public int getTaxaCount() {
        return taxaNum;
    }

    @Override
    public String getTaxaName(int taxaIndex) {
        return taxaNameIndexBiMap.inverse().get(taxaIndex);
    }

    @Override
    public String[] getTaxaNames() {
        String[] taxaNames=new String[taxaList.size()];
        for (int i = 0; i < taxaList.numberOfTaxa(); i++) {
            taxaNames[i]=taxaList.taxaName(i);
        }
        return taxaNames;
    }
}
