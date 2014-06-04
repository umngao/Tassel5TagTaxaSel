/*
 *  HDF5Byte3D
 */
package net.maizegenetics.dna.snp.byte3d;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5Byte3D extends AbstractByte3D {

    private static final int MAX_CACHE_SIZE = 1 << 16;
    private static final int HDF5_BLOCK = 1 << 16;
    private final Map<Long, byte[][]> myDepthCache = new LinkedHashMap<Long, byte[][]>((3 * MAX_CACHE_SIZE) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry<Long, byte[][]> eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };

    private final IHDF5Reader myReader;
    private final int myNumSites;
    private final TaxaList myTaxa;
    private final int myNumAlleles;

    HDF5Byte3D(IHDF5Reader reader) {
        super(6, reader.getIntAttribute(Tassel5HDF5Constants.GENOTYPES_MODULE, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA),
                reader.getIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES)
        );
        myReader = reader;
        myNumSites = reader.getIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        myTaxa = new TaxaListBuilder().buildFromHDF5(reader);
        // TODO - This needs to be retrieved from the HDF5 file
        myNumAlleles = 6;
    }

    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_BLOCK);
    }

    private byte[][] cacheValuesBlock(int taxon, int site, long key) {
        int start = (site / MAX_CACHE_SIZE) * MAX_CACHE_SIZE;
        int realSiteCache = (myNumSites - start < MAX_CACHE_SIZE) ? myNumSites - start : MAX_CACHE_SIZE;
        byte[][] data = myReader.readByteMatrixBlockWithOffset(Tassel5HDF5Constants.getGenotypesDepthPath(myTaxa.taxaName(taxon)), myNumAlleles, realSiteCache, 0, start);
        if (data == null) {
            return null;
        }
        myDepthCache.put(key, data);
        return data;
    }

    @Override
    public byte valueForAllele(int taxon, int site, int allele) {
        long key = getCacheKey(taxon, site);
        byte[][] data = myDepthCache.get(key);
        if (data == null) {
            data = cacheValuesBlock(taxon, site, key);
        }
        return data[allele][site % MAX_CACHE_SIZE];
    }

}
