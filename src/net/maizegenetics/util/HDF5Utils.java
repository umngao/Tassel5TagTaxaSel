package net.maizegenetics.util;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.tag.Tags;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

import java.util.*;

/**
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public final class HDF5Utils {

    private HDF5Utils() {
    }

    public static boolean isTASSEL4HDF5Format(IHDF5Reader reader) {
        return reader.exists(net.maizegenetics.dna.snp.HapMapHDF5Constants.LOCI);
    }

 //TAXA Module
     public static void createHDF5TaxaModule(IHDF5Writer h5w) {
         h5w.createGroup(Tassel5HDF5Constants.TAXA_MODULE);
         h5w.setBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED,false);
     }

    public static void lockHDF5TaxaModule(IHDF5Writer h5w) {
        h5w.setBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED,true);
    }
    
    public static void unlockHDF5TaxaModule(IHDF5Writer h5w) {
        h5w.setBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED,false);
    }

    public static boolean doesTaxaModuleExist(IHDF5Reader reader) {
        return reader.exists(Tassel5HDF5Constants.TAXA_MODULE);
    }

    public static boolean doTaxonCallsExist(IHDF5Reader reader, String taxonName){
        return reader.exists(Tassel5HDF5Constants.getGenotypesCallsPath(taxonName));
    }

    public static boolean doTaxonCallsExist(IHDF5Reader reader, Taxon taxon){
        return reader.exists(Tassel5HDF5Constants.getGenotypesCallsPath(taxon.getName()));
    }

    public static boolean isTaxaLocked(IHDF5Reader reader){
        return reader.getBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED);
    }

    /**
     * Adds a taxon to the taxon module
     * @param h5w
     * @param taxon
     * @return true if new add, false if already exists
     */
    public static boolean addTaxon(IHDF5Writer h5w, Taxon taxon) {
        if(isTaxaLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String path=Tassel5HDF5Constants.getTaxonPath(taxon.getName());
        if(h5w.exists(path)) return false;
        h5w.createGroup(path);
        SetMultimap<String, String> annoMap=taxon.getAnnotationAsMap();
        for (String keys : annoMap.keys()) {
            String s=Joiner.on(",").join(annoMap.get(keys));
            h5w.setStringAttribute(path,keys,s);
        }
        return true;
    }
    
    public static void replaceTaxonAnnotations(IHDF5Writer h5w, Taxon modifiedTaxon) {
        if(isTaxaLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String path=Tassel5HDF5Constants.getTaxonPath(modifiedTaxon.getName());
        if(!h5w.exists(path))  throw new IllegalStateException("Taxon does not already exist to replace annotations of");;
        SetMultimap<String, String> annoMap=modifiedTaxon.getAnnotationAsMap();
        for (String keys : annoMap.keys()) {
            String s=Joiner.on(",").join(annoMap.get(keys));
            h5w.setStringAttribute(path,keys,s);
        }
    }

    public static Taxon getTaxon(IHDF5Reader reader, String taxonName) {
        String taxonPath=Tassel5HDF5Constants.getTaxonPath(taxonName);
        if(!reader.exists(taxonPath)) return null;
        Taxon.Builder tb=new Taxon.Builder(taxonName);
        for (String a : reader.getAllAttributeNames(taxonPath)) {
            for(String s: Splitter.on(",").split(reader.getStringAttribute(taxonPath,a))) {
                tb.addAnno(a,s);
            }
        }
        return tb.build();
    }

    public static List<String> getAllTaxaNames(IHDF5Reader reader) {
        List<String> taxaNames=new ArrayList<>();
        List<HDF5LinkInformation> fields = reader.getAllGroupMemberInformation(Tassel5HDF5Constants.TAXA_MODULE, true);
        for (HDF5LinkInformation is : fields) {
            taxaNames.add(is.getName());
        }
        return taxaNames;
    }


    public static void writeHDF5TaxaNumTaxa(IHDF5Writer h5w, int numTaxa) {
        h5w.setIntAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_NUM_TAXA, numTaxa);
    }

    public static int getHDF5TaxaNumTaxa(IHDF5Reader reader) {
        return reader.getIntAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_NUM_TAXA);
    }


 //GENOTYPE Module

    public static int getHDF5GenotypeTaxaNumber(IHDF5Reader reader){
        return reader.getIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_NUM_TAXA);
    }

    public static void createHDF5GenotypeModule(IHDF5Writer h5w) {
        if(h5w.exists(Tassel5HDF5Constants.GENOTYPES_MODULE)) throw new UnsupportedOperationException("Genotypes module already exists in HDF5 file");
        h5w.createGroup(Tassel5HDF5Constants.GENOTYPES_MODULE);
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED,false);
    }

    public static void lockHDF5GenotypeModule(IHDF5Writer h5w) {
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED,true);
    }
    
    public static void unlockHDF5GenotypeModule(IHDF5Writer h5w) {
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED,false);
    }

    public static boolean doesGenotypeModuleExist(IHDF5Reader reader){
        return reader.exists(Tassel5HDF5Constants.GENOTYPES_MODULE);
    }

    public static boolean isHDF5GenotypeLocked(IHDF5Reader reader){
        if(reader.exists(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH+"/"+Tassel5HDF5Constants.GENOTYPES_LOCKED)==false) return false;
        return reader.getBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED);
    }

    public static void writeHDF5GenotypesMaxNumAlleles(IHDF5Writer h5w, int maxNumAlleles) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, maxNumAlleles);
    }

    public static void writeHDF5GenotypesRetainRareAlleles(IHDF5Writer h5w, boolean retainRareAlleles) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_RETAIN_RARE_ALLELES, retainRareAlleles);
    }

    public static void writeHDF5GenotypesNumTaxa(IHDF5Writer h5w, int numTaxa) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA, numTaxa);
    }

    public static void writeHDF5GenotypesScoreType(IHDF5Writer h5w, String scoreType) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setStringAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, scoreType);
    }

    public static void writeHDF5GenotypesAlleleStates(IHDF5Writer h5w, String[][] aEncodings) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        int numEncodings = aEncodings.length;
        int numStates = aEncodings[0].length;
        MDArray<String> alleleEncodings = new MDArray<>(String.class, new int[]{numEncodings, numStates});
        for (int s = 0; s < numEncodings; s++) {
            for (int x = 0; x < numStates; x++) {
                alleleEncodings.set(aEncodings[s][x], s, x);
            }
        }
        h5w.createStringMDArray(Tassel5HDF5Constants.GENOTYPES_ALLELE_STATES, 100, new int[]{numEncodings, numStates});
        h5w.writeStringMDArray(Tassel5HDF5Constants.GENOTYPES_ALLELE_STATES, alleleEncodings);
    }

    public static byte[] getHDF5GenotypesCalls(IHDF5Reader reader, String taxon) {
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        return reader.readAsByteArray(callsPath);
    }

    public static void writeHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if(h5w.exists(callsPath)) throw new IllegalStateException("Taxa Calls Already Exists:"+taxon);
        h5w.createByteArray(callsPath, calls.length, Math.min(Tassel5HDF5Constants.BLOCK_SIZE,calls.length), Tassel5HDF5Constants.intDeflation);
        writeHDF5EntireArray(callsPath, h5w, calls.length, Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static void replaceHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if(!h5w.exists(callsPath)) throw new IllegalStateException("Taxa Calls Do Not Already Exists to replace");
        writeHDF5EntireArray(callsPath, h5w, calls.length, Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static void replaceHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, int startSite, byte[] calls) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if(!h5w.exists(callsPath)) throw new IllegalStateException("Taxa Calls Do Not Already Exists to replace");
        if(startSite%Tassel5HDF5Constants.BLOCK_SIZE!=0) throw new IllegalStateException("Taxa Calls Start Site not a multiple of the block size");
        writeHDF5Block(callsPath,h5w,Tassel5HDF5Constants.BLOCK_SIZE,startSite/Tassel5HDF5Constants.BLOCK_SIZE,calls);
    }

    public static byte[][] getHDF5GenotypesDepth(IHDF5Reader reader, String taxon) {
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if(reader.exists(callsPath)) {
            return reader.readByteMatrix(callsPath);
        } else {
            return null;
        }
    }

    public static boolean doesGenotypeDepthExist(IHDF5Reader reader) {
        List<String> taxaNames=getAllTaxaNames(reader);
        boolean depthExist=false;
        for (String taxon : taxaNames) {
            String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
            if(reader.exists(callsPath)) {
                depthExist=true;
                break;
            }
        }
        return depthExist;
    }

    public static void writeHDF5GenotypesDepth(IHDF5Writer h5w, String taxon, byte[][] depth) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if(h5w.exists(callsPath)) throw new IllegalStateException("Taxa Depth Already Exists:"+taxon);
        h5w.createByteMatrix(callsPath, depth.length, depth[0].length, 6, Math.min(Tassel5HDF5Constants.BLOCK_SIZE,depth[0].length), Tassel5HDF5Constants.intDeflation);
        writeHDF5EntireArray(callsPath, h5w, depth[0].length, Tassel5HDF5Constants.BLOCK_SIZE, depth);
    }

    public static void replaceHDF5GenotypesDepth(IHDF5Writer h5w, String taxon, byte[][] depth) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if(!h5w.exists(callsPath)) throw new IllegalStateException("Taxa Depth Does Not Already Exists to Replace");
        writeHDF5EntireArray(callsPath, h5w, depth[0].length, Tassel5HDF5Constants.BLOCK_SIZE, depth);
    }

    public static byte[] getHDF5Alleles(IHDF5Reader reader, WHICH_ALLELE allele) {
        return reader.readByteMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 1, getHDF5PositionNumber(reader),
                (long)allele.index(), 0)[0];
    }




//    Positions/numSites
    public static int getHDF5PositionNumber(IHDF5Reader reader){
        return reader.getIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH,Tassel5HDF5Constants.POSITION_NUM_SITES);
    }

    public static void createHDF5PositionModule(IHDF5Writer h5w) {
        h5w.createGroup(Tassel5HDF5Constants.POSITION_MODULE);
    }

    public static void writeHDF5PositionNumSite(IHDF5Writer h5w, int numSites) {
        h5w.setIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES, numSites);
    }

    public static byte[] getHDF5ReferenceAlleles(IHDF5Reader reader) {
        return getHDF5Alleles(reader,Tassel5HDF5Constants.REF_ALLELES, 0, getHDF5PositionNumber(reader));
    }

    public static byte[] getHDF5ReferenceAlleles(IHDF5Reader reader, int startSite, int numSites) {
        return getHDF5Alleles(reader,Tassel5HDF5Constants.REF_ALLELES, startSite, numSites);
    }

    public static byte[] getHDF5AncestralAlleles(IHDF5Reader reader) {
        return getHDF5Alleles(reader,Tassel5HDF5Constants.ANC_ALLELES, 0, getHDF5PositionNumber(reader));
    }

    public static byte[] getHDF5AncestralAlleles(IHDF5Reader reader, int startSite, int numSites) {
        return getHDF5Alleles(reader,Tassel5HDF5Constants.ANC_ALLELES, startSite, numSites);
    }


    private static byte[] getHDF5Alleles(IHDF5Reader reader, String allelePath, int startSite, int numSites) {
        if(reader.exists(allelePath)) {
            return reader.readByteArrayBlockWithOffset(allelePath, numSites, startSite);
        }
        byte[] unknown=new byte[numSites];
        Arrays.fill(unknown, GenotypeTable.UNKNOWN_ALLELE);
        return unknown;
    }


    //Tags Module
    public static void createHDF5TagModule(IHDF5Writer h5w, int tagLengthInLong) {
        if(h5w.exists(Tassel5HDF5Constants.TAG_MODULE)) throw new UnsupportedOperationException("Tag module already exists in HDF5 file");
        h5w.createGroup(Tassel5HDF5Constants.TAG_MODULE);
        h5w.setBooleanAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAG_LOCKED,false);
        h5w.setIntAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LENGTH_LONG, tagLengthInLong);
        h5w.setIntAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_COUNT, 0);
    }

    public static boolean isHDF5TagLocked(IHDF5Reader reader) {
        if (reader.exists(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH + "/" + Tassel5HDF5Constants.TAG_LOCKED) == false)
            return false;
        return reader.getBooleanAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LOCKED);
    }

    public static int getHDF5TagCount(IHDF5Reader reader){
        return reader.getIntAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_COUNT);
    }

    public static int getHDF5TagLengthInLong(IHDF5Reader reader){
        return reader.getIntAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LENGTH_LONG);
    }

    public static boolean doTagsExist(IHDF5Reader reader){
        return reader.exists(Tassel5HDF5Constants.TAGS);
    }

    public static void writeHDF5Tags(IHDF5Writer h5w, Tags tags) {
        if(doTagsExist(h5w)) throw new UnsupportedOperationException("Tags already exists in HDF5 file");
        if(isHDF5TagLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        int blockSize = Math.min(tags.getTagCount(), Tassel5HDF5Constants.BLOCK_SIZE);
        h5w.createLongMatrix(Tassel5HDF5Constants.TAGS, tags.getTagSizeInLong(), tags.getTagCount(), tags.getTagSizeInLong(),
                blockSize, Tassel5HDF5Constants.intDeflation);
        h5w.createIntArray(Tassel5HDF5Constants.TAG_LENGTHS, tags.getTagCount(), blockSize,Tassel5HDF5Constants.intDeflation);
        h5w.setIntAttribute(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_COUNT, tags.getTagCount());
        int blocks = ((tags.getTagCount() - 1) / blockSize) + 1;
        for (int block = 0; block < blocks; block++) {
            int startPos = block * blockSize;
            int length = Math.min(tags.getTagCount() - startPos, blockSize);
            long[][] sval = new long[tags.getTagSizeInLong()][length];
            int[] tLenBlock = new int[length];
            for (int k = 0; k < length; k++) {
                long[] tg=tags.getTag(startPos+k);
                for (int j = 0; j < tg.length; j++) {
                    sval[j][k]=tg[j];
                }
                tLenBlock[k]=tags.getTagLength(k);
            }
            writeHDF5Block(Tassel5HDF5Constants.TAGS, h5w, blockSize, block, sval);
            writeHDF5Block(Tassel5HDF5Constants.TAG_LENGTHS, h5w, blockSize, block, tLenBlock);
        }
    }

    public static long[][] getTags(IHDF5Reader reader) {
        return reader.readLongMatrix(Tassel5HDF5Constants.TAGS);
    }

    public static int[] getTagLengths(IHDF5Reader reader) {
        return reader.readIntArray(Tassel5HDF5Constants.TAG_LENGTHS);
    }

    public static void createHDF5TagByTaxaDist(IHDF5Writer h5w, boolean isTaxaDirection, TaxaList taxaList) {
        if(h5w.exists(Tassel5HDF5Constants.TAG_DIST)) throw new UnsupportedOperationException("TagDist module already exists in HDF5 file");
        if(isTaxaDirection) {
            h5w.createByteMatrix(Tassel5HDF5Constants.TAG_DIST, 0, getHDF5TagCount(h5w),
                    1, Math.min(getHDF5TagCount(h5w), Tassel5HDF5Constants.BLOCK_SIZE), Tassel5HDF5Constants.intDeflation);
            h5w.setBooleanAttribute(Tassel5HDF5Constants.TAG_DIST,Tassel5HDF5Constants.TAG_DIST_CHUNK,true);
        } else {
            h5w.createByteMatrix(Tassel5HDF5Constants.TAG_DIST, taxaList.size(), getHDF5TagCount(h5w),
                    taxaList.size(), 1, Tassel5HDF5Constants.intDeflation);
            h5w.setBooleanAttribute(Tassel5HDF5Constants.TAG_DIST,Tassel5HDF5Constants.TAG_DIST_CHUNK,false);
            for (int t = 0; t < taxaList.size(); t++) {
                h5w.setStringAttribute(Tassel5HDF5Constants.TAG_DIST,"TN"+t,taxaList.taxaName(t));
            }
        }
    }

    public static int addTaxonTagDistribution(IHDF5Writer h5w, Taxon taxon, byte[] dist) {
        if(!isTagsByTaxaInTaxaDirection(h5w)) throw new IllegalStateException("Chunking in wrong direction for adding taxa");
        addTaxon(h5w,taxon);
        int tbtTaxaNum=(int)h5w.getDataSetInformation(Tassel5HDF5Constants.TAG_DIST).getDimensions()[0];
        //System.out.println(h5w.getDataSetInformation(Tassel5HDF5Constants.TAG_DIST).toString());
        h5w.writeByteMatrixBlock(Tassel5HDF5Constants.TAG_DIST, new byte[][]{dist}, tbtTaxaNum, 0);
        h5w.setStringAttribute(Tassel5HDF5Constants.TAG_DIST,"TN"+tbtTaxaNum,taxon.getName());
        return tbtTaxaNum;
    }

    public static byte[] getTagDistForTaxon(IHDF5Reader reader, int taxonIndex) {
        if(!isTagsByTaxaInTaxaDirection(reader)) throw new IllegalStateException("Chunking in wrong direction for reading taxa");
        //getting read count and checking read direction may slow things.
        return reader.readByteMatrixBlock(Tassel5HDF5Constants.TAG_DIST,1,getHDF5TagCount(reader),taxonIndex,0l)[0];
    }

    public static boolean doTagsByTaxaExist(IHDF5Reader reader){
        return reader.exists(Tassel5HDF5Constants.TAG_DIST);
    }

    public static boolean isTagsByTaxaInTaxaDirection(IHDF5Reader reader){
        return reader.getBooleanAttribute(Tassel5HDF5Constants.TAG_DIST, Tassel5HDF5Constants.TAG_DIST_CHUNK);
    }

    public static int getNumberOfTaxaInTBT(IHDF5Reader reader){
        return (int)reader.getDataSetInformation(Tassel5HDF5Constants.TAG_DIST).getDimensions()[0];
    }

    public static Map<String, Integer> getTBTMapOfRowIndices(IHDF5Reader reader) {
        int tbtTaxaNum=(int)reader.getDataSetInformation(Tassel5HDF5Constants.TAG_DIST).getDimensions()[0];
        Map<String, Integer> taxaNameToRowIndex=new HashMap<>(tbtTaxaNum);
        for (int t = 0; t < tbtTaxaNum; t++) {
            taxaNameToRowIndex.put(reader.getStringAttribute(Tassel5HDF5Constants.TAG_DIST,"TN"+t),t);
        }
        return taxaNameToRowIndex;
    }

//    public static final String TAG_COUNT = "tagCount";
//    public static final String TAG_LENGTH_LONG = "tagLengthLong";
//    public static final String TAG_LOCKED = "locked";
//    public static final String TAGS = TAG_MODULE + "/Tags";
//    public static final String TAG_LENGTHS = TAG_MODULE + "/TagLength";
//    public static final String TAG_DIST = TAG_MODULE + "/TagDist";
//    public static final String TAG_DIST_CHUNK = "chunkDirection";

  //Writers for these should also be implemented, but there are some data scale issue that they are written in blocks.
    //see public PositionListBuilder(IHDF5Writer h5w, PositionList a)
    //    Positions/Positions
//    Positions/Chromosome
//    Positions/ChromosomeIndices
//    Positions/SnpIds



    /**
     *
     * @param objectPath
     * @param myWriter
     * @param objMaxLength
     * @param blockSize
     * @param val
     */
    public static void writeHDF5EntireArray(String objectPath, IHDF5Writer myWriter, int objMaxLength, int blockSize, Object val) {
        int blocks = ((objMaxLength - 1) / blockSize) + 1;
        for (int block = 0; block < blocks; block++) {
            int startPos = block * blockSize;
            int length = Math.min(objMaxLength - startPos, blockSize);
            if (val instanceof byte[][]) {
                byte[][] oval = (byte[][]) val;
                byte[][] sval = new byte[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof int[][]) {
                int[][] oval = (int[][]) val;
                int[][] sval = new int[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof long[][]) {
                long[][] oval = (long[][]) val;
                long[][] sval = new long[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof byte[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((byte[]) val, startPos, startPos + length));
            } else if (val instanceof float[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((float[]) val, startPos, startPos + length));
            } else if (val instanceof int[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((int[]) val, startPos, startPos + length));
            } else if (val instanceof String[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((String[]) val, startPos, startPos + length));
            }
        }
    }

    /**
     *
     * @param objectPath
     * @param myWriter
     * @param blockSize
     * @param block
     * @param val
     */
    public static void writeHDF5Block(String objectPath, IHDF5Writer myWriter, int blockSize, int block, Object val) {
        int startPos = block * blockSize;
        if (val instanceof byte[][]) {
            byte[][] bval = (byte[][]) val;
            myWriter.writeByteMatrixBlockWithOffset(objectPath, bval, bval.length, bval[0].length, 0l, (long) startPos);
        } else if (val instanceof byte[]) {
            byte[] fval = (byte[]) val;
            myWriter.writeByteArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof float[]) {
            float[] fval = (float[]) val;
            myWriter.writeFloatArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[]) {
            int[] fval = (int[]) val;
            myWriter.writeIntArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[][]) {
            int[][] ival = (int[][]) val;
            myWriter.writeIntMatrixBlockWithOffset(objectPath, ival, ival.length, ival[0].length, 0l, (long) startPos);
        } else if (val instanceof long[][]) {
            long[][] lval = (long[][]) val;
            myWriter.writeLongMatrixBlockWithOffset(objectPath, lval, lval.length, lval[0].length, 0l, (long) startPos);
        } else  if (val instanceof String[]) {
            String[] sval = (String[]) val;
            myWriter.writeStringArrayBlockWithOffset(objectPath, sval, sval.length, (long) startPos);
        }
    }



}
