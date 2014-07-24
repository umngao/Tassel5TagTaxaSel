/*
 *  HDF5ByteGenotype
 */
package net.maizegenetics.dna.snp.genotypecall;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.concurrent.ExecutionException;

/**
 * HDF5 implementation of GenotypeTable. Uses caching of GenotypeTable,
 * alleleCounts, MAF, and siteCoverage
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
class HDF5ByteGenotypeCallTable extends AbstractGenotypeCallTable {

    private static final int SHIFT_AMOUNT = 16;

    private final String[] genotypePaths;
    /**
     * Byte representations of DNA sequences are stored in blocks of 65536 sites
     */
    private static final int HDF5_GENOTYPE_BLOCK_SIZE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(HDF5_GENOTYPE_BLOCK_SIZE - 1);
    private final IHDF5Reader myHDF5Reader;

    private final LoadingCache<Long, byte[]> myGenoCache;
    private final CacheLoader<Long, byte[]> myGenoLoader = new CacheLoader<Long, byte[]>() {
        @Override
        public byte[] load(Long key) {
            long offset = getSiteStartFromKey(key) << SHIFT_AMOUNT;
            byte[] data;
            synchronized (myHDF5Reader) {
                data = myHDF5Reader.readAsByteArrayBlockWithOffset(getTaxaGenoPath(getTaxonFromKey(key)), HDF5_GENOTYPE_BLOCK_SIZE, offset);
            }
            return data;
        }
    };

    private final LoadingCache<Integer, SiteBlockAttr> mySiteAnnoCache; //key = site
    private final CacheLoader<Integer, SiteBlockAttr> siteAnnotLoader = new CacheLoader<Integer, SiteBlockAttr>() {
        int lastCachedStartSite = Integer.MIN_VALUE;
        int[][] af;
        byte[][] afOrder;
        float[] maf;
        float[] paf;

        @Override
        public SiteBlockAttr load(Integer key) {
            int startSite = getStartSite(key);
            int length = Math.min(HDF5_GENOTYPE_BLOCK_SIZE, numberOfSites() - startSite);
            System.out.println("Reading from HDF5 site anno:" + startSite);
            System.out.println("");
            synchronized (myHDF5Reader) {
                af = myHDF5Reader.readIntMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_CNT, 6, length, 0l, startSite);
                afOrder = myHDF5Reader.readByteMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 6, length, 0l, startSite);
                maf = myHDF5Reader.readFloatArrayBlockWithOffset(Tassel5HDF5Constants.MAF, length, startSite);
                paf = myHDF5Reader.readFloatArrayBlockWithOffset(Tassel5HDF5Constants.SITECOV, length, startSite);
                lastCachedStartSite = startSite;
            }
            return new SiteBlockAttr(startSite, afOrder, af, maf, paf);
        }
    };

    private class SiteBlockAttr {

        private final int startSite;  //4
        private final byte[][] myAlleleFreqOrder;  //[
        private final int[][] myAlleleCnt;  //2-6*4=24,  sorted by allele frequency of myAlleleFreqOrder
        private final float[] maf; //4
        private final float[] siteCov;  //4

        public SiteBlockAttr(int startSite, byte[][] myAlleleFreqOrder, int[][] myAlleleCnt,
                float[] maf, float[] siteCov) {
            this.startSite = startSite;
            this.myAlleleFreqOrder = myAlleleFreqOrder;
            this.myAlleleCnt = myAlleleCnt;
            this.maf = maf;
            this.siteCov = siteCov;
        }

        public int[][] getAllelesSortedByFrequency(int site) {
            int offset = site - startSite;
            int alleleCnt = 0;
            while (myAlleleFreqOrder[alleleCnt][offset] != GenotypeTable.UNKNOWN_ALLELE) {
                alleleCnt++;
            }
            int result[][] = new int[2][alleleCnt];
            for (int i = 0; i < alleleCnt; i++) {
                result[0][i] = myAlleleFreqOrder[i][offset];
                result[1][i] = myAlleleCnt[result[0][i]][offset];
            }
            return result;
        }

        public float getMAF(int site) {
            return maf[site - startSite];
        }

        public float getSiteCoverage(int site) {
            return siteCov[site - startSite];
        }

    }

    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_GENOTYPE_BLOCK_SIZE);
    }

    private static int getTaxonFromKey(long key) {
        return (int) (key >>> 33);
    }

    private static int getSiteStartFromKey(long key) {
        return (int) ((key << 33) >>> 33);
    }

    private static int getStartSite(int site) {
        return site & SITE_BLOCK_MASK;
    }

    private String getTaxaGenoPath(int taxon) {
        return genotypePaths[taxon];
    }

    private HDF5ByteGenotypeCallTable(IHDF5Reader reader, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        super(numTaxa, numSites, phased, alleleEncodings);
        genotypePaths = new String[numTaxa];
        TaxaList tL = new TaxaListBuilder().buildFromHDF5Genotypes(reader);  //not the most efficient thing to do, but ensures sort is the same.
        for (int i = 0; i < numTaxa; i++) {
            genotypePaths[i] = Tassel5HDF5Constants.getGenotypesCallsPath(tL.taxaName(i));
        }
        myHDF5Reader = reader;
        myGenoCache = CacheBuilder.newBuilder()
                .maximumSize(Runtime.getRuntime().availableProcessors() * 3)
                .build(myGenoLoader);
        mySiteAnnoCache = CacheBuilder.newBuilder()
                .maximumSize(150)
                .build(siteAnnotLoader);
    }

    static HDF5ByteGenotypeCallTable getInstance(IHDF5Reader reader) {
        int numTaxa = reader.getIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA);
        int numSites = reader.getIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        String[][] alleleEncodings = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
        return new HDF5ByteGenotypeCallTable(reader, numTaxa, numSites, false, alleleEncodings);
    }

    @Override
    public byte genotype(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        try {
            byte[] data = myGenoCache.get(key);
            return data[site % HDF5_GENOTYPE_BLOCK_SIZE];
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new IllegalStateException("HDF5ByteGenotype: getBase: Error getting base from cache.");
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype(taxon, site));
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        try {
            SiteBlockAttr sa = mySiteAnnoCache.get(getStartSite(site));
            return sa.getAllelesSortedByFrequency(site);
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new UnsupportedOperationException("Error in getMinorAlleleFrequency from cache");
        }
    }

    @Override
    public double minorAlleleFrequency(int site) {
        try {
            SiteBlockAttr sa = mySiteAnnoCache.get(getStartSite(site));
            return sa.getMAF(site);
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new UnsupportedOperationException("Error in getMinorAlleleFrequency from cache");
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
