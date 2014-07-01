/*
 * TagsByTaxaByteHDF5TaxaGroups
 */
package net.maizegenetics.dna.tag;

import cern.colt.list.IntArrayList;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import com.google.common.collect.BiMap;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;

import static net.maizegenetics.util.Tassel5HDF5Constants.*;

/**
 * Tags by Taxa (TBT) Builder based on the HDF5 data structure.  This version uses the standard 2-dimension matrix with
 * chunking either set for taxa or tags.  Normally, the data structure will be built by adding taxa to the structure.
 * <p></p>
 * Discovery of SNPs requires that the TBT is pivoted to the tags based chunking, which is done by calling the createTagOriented.
 * 
 * @author  Ed Buckler
 */
public class TagsByTaxaHDF5Builder {

    //private final Chunking cnkDir;

//    private int tagCount = 0;
//    private final IHDF5Writer h5;
//    private final BiMap<String, Integer> taxonNameToIndexMap;

//    /**
//     * Create Tags by taxa module within a HDF5 file
//     * @param newHDF5file file name
//     * @param tags the list of tags to be use
//     * @return
//     */
//    public static final TagsByTaxaHDF5Builder createTaxaIncremental(String newHDF5file, Tags tags) {
//        return new TagsByTaxaHDF5Builder(newHDF5file, Chunking.Taxa,tags, null);
//    }

//    /**
//     * Used to add taxa to an existing Tags by Taxa module
//     * @param existingHDF5file
//     * @return
//     */
//    public static final TagsByTaxaHDF5Builder openTaxaIncremental(String existingHDF5file) {
//        return new TagsByTaxaHDF5Builder(existingHDF5file,Chunking.Taxa,null,null);
//    }
//
//    public static final TagsByTaxaHDF5Builder createTagOriented(String newHDF5file, Tags tags, TaxaList taxaList) {
//        return new TagsByTaxaHDF5Builder(newHDF5file,Chunking.Tags,tags,taxaList);
//    }

    public static final void create(String newHDF5file, Map<Tag,TaxaDistribution> tagTaxaMap, TaxaList taxaList) {
        IHDF5Writer h5w = HDF5Factory.configure(new File(newHDF5file))
                .useUTF8CharacterEncoding().writer();
        Map.Entry<Tag, TaxaDistribution> firstEntry=tagTaxaMap.entrySet().iterator().next();
        Tag aTag=firstEntry.getKey();
        int maxTaxa=firstEntry.getValue().maxTaxa();
        if(maxTaxa!=taxaList.numberOfTaxa()) throw new IllegalStateException("Taxa Distribution does not agree with size of TaxaList");
        TaxaListBuilder.createHDF5TaxaList(h5w,taxaList);
        HDF5Utils.createHDF5TagModule(h5w,aTag.seq2Bit().length);
        ArrayList<Map.Entry<Tag, TaxaDistribution>>[] blockList=new ArrayList[Tassel5HDF5Constants.TAGS_BIN_NUM];
        for (int i = 0; i < blockList.length; i++) {
            blockList[i]=new ArrayList<>();
        }
        for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaMap.entrySet()) {
            blockList[entry.getKey().hashCode()>>>Tassel5HDF5Constants.HASH_SHIFT_TO_TAG_BIN].add(entry);
        }
        //todo parallelize this
        int numThreads=Runtime.getRuntime().availableProcessors();
        ExecutorService pool= Executors.newFixedThreadPool(numThreads - 1);
        List<WriteTagBucket> bucketList=new ArrayList<>();
        for (int i = 0; i < blockList.length; i++) {
            //System.out.println("Using File:"+inputSeqFile.toString());
            bucketList.add(new WriteTagBucket(h5w,i,blockList[i]));
        }
        try{pool.invokeAll(bucketList);}
        catch (Exception ie) {
            ie.printStackTrace();
        }
    }

    public static final TagsByTaxa openForReading(String aHDF5file) {
        throw new UnsupportedOperationException("Not implemented yet.");
    }

}

class WriteTagBucket implements Callable<Integer> {
    private final int bucket;
    private final ArrayList<Map.Entry<Tag, TaxaDistribution>> tagDistList;
    private final IHDF5Writer h5w;
    private final int maxTaxa=192;

    WriteTagBucket(IHDF5Writer h5w, int bucket, ArrayList<Map.Entry<Tag, TaxaDistribution>> tagDistList) {
        this.bucket=bucket;
        this.tagDistList=tagDistList;
        this.h5w=h5w;
    }

    @Override
    public Integer call() throws Exception {
        if(tagDistList.size()<1) return bucket;
        int block=tagDistList.get(0).getKey().hashCode()>>>Tassel5HDF5Constants.HASH_SHIFT_TO_TAG_BIN;
        long[][] tags=new long[tagDistList.get(0).getKey().seq2Bit().length][tagDistList.size()];
        short[] length=new short[tagDistList.size()];
        IntArrayList taxaDist=new IntArrayList(tagDistList.size());
        IntArrayList tagDistOffset=new IntArrayList(tagDistList.size());
        int t_i=0;
        for (Map.Entry<Tag, TaxaDistribution> entry : tagDistList) {
            Tag t=entry.getKey();
            int li=0; for (long l : t.seq2Bit()) {tags[li++][t_i]=l;}
            length[t_i]=t.seqLength();
            int[] encodeTagDist=entry.getValue().encodeTaxaDepth();
            tagDistOffset.add(taxaDist.size());
            for (int i : encodeTagDist) {
                taxaDist.add(i);
            }
            t_i++;
        }
        tagDistOffset.add(taxaDist.size());
        taxaDist.trimToSize();
        tagDistOffset.trimToSize();
        //TODO use snappy for compression over gzip for speed and greater parallelization in java
//        byte[] compressed = Snappy.compress(input.getBytes("UTF-8"));
//        byte[] uncompressed = Snappy.uncompress(compressed);
        HDF5Utils.writeTagDistributionBucket(h5w,block, tags,length,taxaDist.elements(),maxTaxa,tagDistOffset.elements());
        return bucket;
    }
}
