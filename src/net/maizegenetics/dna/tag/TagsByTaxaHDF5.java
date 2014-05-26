/*
 * TagsByTaxaByteHDF5TaxaGroups
 */
package net.maizegenetics.dna.tag;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.util.*;

import java.util.*;
import java.util.concurrent.ExecutionException;

/**
 * Tags by Taxa file based on the HDF5 data structure.  This version is optimized
 * for rapid access of tags within taxa (ie it buffers the tag counts within one
 * taxon).  It is good for adding, removing, and combining taxa
 * 
 * @author edbuckler
 */
public class TagsByTaxaHDF5 extends AbstractTagsByTaxa {
    public enum Chunking {
        Taxa, Tags
    }
    private boolean inTaxaChunking;

    //Guava cache loader
    //load upto 80000 tags or

    static int chunkSize = 1 << 16;
    int tagCount = 0;
    int tagChunks = 0;
    IHDF5Reader reader = null;
    List<String> taxaDirList, taxaNameList;
    Map<String, String> taxaNameDirTreeMap;

    private LoadingCache<Integer, byte[]> taxaCache=null;

    private CacheLoader<Integer, byte[]> taxaLoader = new CacheLoader<Integer, byte[]>() {
        public byte[] load(Integer key) {
           return HDF5Utils.getTagDistForTaxon(reader,key); //HDF5Utils.
        }
    };




    public TagsByTaxaHDF5(String infile) {
        reader = HDF5Factory.openForReading(infile);
        tagCount = HDF5Utils.getHDF5TagCount(reader);
        inTaxaChunking= HDF5Utils.isTagsByTaxaInTaxaDirection(reader);
        tagLengthInLong = HDF5Utils.getHDF5TagLengthInLong(reader);
        taxaNum = HDF5Utils.getNumberOfTaxaInTBT(reader);
        this.tags = HDF5Utils.getTags(reader);

        //todo need to make tag lengths an integer
        int[] tL=HDF5Utils.getTagLengths(reader);
        this.tagLength=new byte[tL.length];
        for (int i = 0; i < tagLength.length; i++) {
            tagLength[i]=(byte)tL[i];
        }



        if(inTaxaChunking) {
            taxaCache = CacheBuilder.newBuilder()
//                    .maximumWeight(300_000_000)
//                    .weigher(weighByLength)
                    .build(taxaLoader);
        }
    }




    @Override
    public int getIndexOfTaxaName(String taxon) {
        int index = Collections.binarySearch(taxaNameList, taxon);
        return index;
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        return getReadCountDistributionForTaxon(taxaIndex)[tagIndex];
    }

    public byte[] getReadCountDistributionForTaxon(int taxaIndex) {
        try {
            return taxaCache.get(taxaIndex);
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
        return taxaNameList.get(taxaIndex);
    }

    @Override
    public String[] getTaxaNames() {
        String[] array = taxaNameList.toArray(new String[taxaNameList.size()]);
        return array;
    }
}
