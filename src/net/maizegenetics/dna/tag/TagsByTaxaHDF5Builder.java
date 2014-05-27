/*
 * TagsByTaxaByteHDF5TaxaGroups
 */
package net.maizegenetics.dna.tag;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import com.google.common.collect.BiMap;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import static net.maizegenetics.dna.tag.TagsByTaxaHDF5.Chunking;

import java.io.File;
import java.util.*;

/**
 * Tags by Taxa (TBT) Builder based on the HDF5 data structure.  This version uses the standard 2-dimension matrix with
 * chunking either set for taxa or tags.  Normally, the data structure will be built by adding taxa to the structure.
 * <p></p>
 * Discovery of SNPs requires that the TBT is pivoted to the tags based chunking, which is done by calling the createTagOriented.
 * 
 * @author  Ed Buckler
 */
public class TagsByTaxaHDF5Builder {

    private final Chunking cnkDir;

    private int tagCount = 0;
    private final IHDF5Writer h5;
    private final BiMap<String, Integer> taxonNameToIndexMap;

    /**
     * Create Tags by taxa module within a HDF5 file
     * @param newHDF5file file name
     * @param tags the list of tags to be use
     * @return
     */
    public static final TagsByTaxaHDF5Builder createTaxaIncremental(String newHDF5file, Tags tags) {
        return new TagsByTaxaHDF5Builder(newHDF5file, Chunking.Taxa,tags, null);
    }

    /**
     * Used to add taxa to an existing Tags by Taxa module
     * @param existingHDF5file
     * @return
     */
    public static final TagsByTaxaHDF5Builder openTaxaIncremental(String existingHDF5file) {
        return new TagsByTaxaHDF5Builder(existingHDF5file,Chunking.Taxa,null,null);
    }

    public static final TagsByTaxaHDF5Builder createTagOriented(String newHDF5file, Tags tags, TaxaList taxaList) {
        return new TagsByTaxaHDF5Builder(newHDF5file,Chunking.Tags,tags,taxaList);
    }

//    public static final TagsByTaxaHDF5Builder createTagOriented(String srcHDF5File, String dstHDF5file) {
//        return new TagsByTaxaHDF5Builder(newHDF5file,BuildDirection.AddTaxa,tags);
//    }

    private TagsByTaxaHDF5Builder(String theHDF5file, TagsByTaxaHDF5.Chunking chunkDirection, Tags tags, TaxaList taxaList) {
        cnkDir=chunkDirection;
        h5 = HDF5Factory.configure(new File(theHDF5file))
                .useUTF8CharacterEncoding().writer();
        if(tags==null) {
            if(!HDF5Utils.doTagsExist(h5)) throw new IllegalStateException("File :"+theHDF5file+" does not have the needed tags" +
                            " to build a TBT file.");
        } else {
            if(HDF5Utils.doesTaxaModuleExist(h5)==false) HDF5Utils.createHDF5TaxaModule(h5);
            HDF5Utils.createHDF5TagModule(h5,tags.getTagSizeInLong());
            HDF5Utils.writeHDF5Tags(h5,tags);
        }
        tagCount=HDF5Utils.getHDF5TagCount(h5);
        if(HDF5Utils.doTagsByTaxaExist(h5)) {
            //validate the taxa list
        } else {
            HDF5Utils.createHDF5TagByTaxaDist(h5,chunkDirection==Chunking.Taxa, taxaList);
        }
        taxonNameToIndexMap=HDF5Utils.getTBTMapOfRowIndices(h5);
    }

    public TagsByTaxaHDF5Builder addTag(long[] tag, byte[] dist) {
        //check if Tags
        if(cnkDir!=Chunking.Tags) throw new IllegalStateException("This builder with these options can only add Tags");
        //do something
        return this;
    }

    public TagsByTaxa build() {
        h5.close();
        return new TagsByTaxaHDF5(h5.getFile().getAbsolutePath());
    }


    public TagsByTaxaHDF5Builder addTaxon(Taxon taxon, byte[] values) {
        if(cnkDir!=Chunking.Taxa) throw new IllegalStateException("This builder with these options can only add Taxa");
        synchronized (h5) {
            if (values.length != this.tagCount) {
                System.err.printf("Taxon (%s) does not have the right number of sites (%d)%n", taxon.getName(), values.length);
                throw new IllegalStateException(String.format("Taxon (%s) does not have the right number of sites (%d)%n", taxon.getName(), values.length));
            }
            if(taxonNameToIndexMap.containsKey(taxon.getName())) {
                throw new IllegalStateException(String.format("Taxon (%s) already exists in the TBT table, please use replace if this is not an error%n", taxon.getName()));
            }
            int rowIndex=HDF5Utils.addTaxonTagDistribution(h5,taxon,values);
            taxonNameToIndexMap.put(taxon.getName(),rowIndex);
            return this;
        }
    }

    public TagsByTaxaHDF5Builder replaceTaxon(Taxon oldTaxon, Taxon newTaxon, byte[] values) {
        synchronized (h5) {
//            if (values.length != this.tagCount) {
//                System.err.printf("Taxon (%s) does not have the right number of sites (%d)%n", taxon.getName(), values.length);
//                throw new IllegalStateException(String.format("Taxon (%s) does not have the right number of sites (%d)%n", taxon.getName(), values.length));
//            }
//            if(taxonNameToIndexMap.containsKey(taxon.getName())) {
//                throw new IllegalStateException(String.format("Taxon (%s) already exists in the TBT table, please use replace if this is not an error%n", taxon.getName()));
//            }
//            int rowIndex=HDF5Utils.addTaxonTagDistribution(h5,taxon,values);
//            taxonNameToIndexMap.put(taxon.getName(),rowIndex);
            return this;
        }
    }

}
