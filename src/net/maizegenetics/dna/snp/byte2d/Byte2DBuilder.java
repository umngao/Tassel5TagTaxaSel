/*
 *  Byte2DBuilder
 */
package net.maizegenetics.dna.snp.byte2d;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 * Builder to store 3 dimensional byte encoded values.
 *
 * @author Terry Casstevens
 */
public class Byte2DBuilder {

    private SuperByteMatrix myValues = null;
    private final boolean myIsHDF5;
    private IHDF5Writer myHDF5Writer = null;
    private final int myNumSites;
    private int myNumTaxa = 0;
    private final GenotypeTable.SITE_SCORE_TYPE mySiteScoreType;

    private Byte2DBuilder(IHDF5Writer writer, int numSites, GenotypeTable.SITE_SCORE_TYPE siteScoreType) {
        myIsHDF5 = true;
        myHDF5Writer = writer;
        myNumSites = numSites;
        mySiteScoreType = siteScoreType;
    }

    private Byte2DBuilder(int numTaxa, int numSites, GenotypeTable.SITE_SCORE_TYPE siteScoreType) {
        myIsHDF5 = false;
        myHDF5Writer = null;
        myNumSites = numSites;
        myNumTaxa = numTaxa;
        mySiteScoreType = siteScoreType;
        myValues = SuperByteMatrixBuilder.getInstance(myNumTaxa, myNumSites);
    }

    public static Byte2DBuilder getInstance(int numTaxa, int numSites, GenotypeTable.SITE_SCORE_TYPE siteScoreType) {
        return new Byte2DBuilder(numTaxa, numSites, siteScoreType);
    }

    public static Byte2DBuilder getInstance(IHDF5Writer writer, int numSites, GenotypeTable.SITE_SCORE_TYPE siteScoreType) {
        return new Byte2DBuilder(writer, numSites, siteScoreType);
    }

    public static FilterByte2D getFilteredInstance(Byte2D base, FilterGenotypeTable filterGenotypeTable) {
        return new FilterByte2D(base, filterGenotypeTable);
    }

    /**
     * Set values for the all sites for a taxon simultaneously.
     *
     * @param taxon Index of taxon
     * @param siteOffset site offset
     * @param values array[sites] of all values
     *
     * @return builder
     */
    public Byte2DBuilder setRangeForTaxon(int taxon, int siteOffset, byte[] values) {
        if (myIsHDF5) {
            throw new IllegalStateException("Byte3DBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        for (int s = 0; s < values.length; s++) {
            set(taxon, s + siteOffset, values[s]);
        }
        return this;
    }

    /**
     * Set value for taxon and site. Value should have already been translated.
     *
     * @param taxon taxon
     * @param site site
     * @param value value
     *
     * @return builder
     */
    public Byte2DBuilder set(int taxon, int site, byte value) {
        if (myIsHDF5) {
            throw new IllegalStateException("Byte3DBuilder: set: use addTaxon for HDF5 files.");
        }
        myValues.set(taxon, site, value);
        return this;
    }

    /**
     * Set values for the all sites for a taxon simultaneously. Values should
     * have already been translated.
     *
     * @param taxon Index of taxon
     * @param values array[sites] of all values
     *
     * @return builder
     */
    public Byte2DBuilder set(int taxon, byte[] values) {
        if (myIsHDF5) {
            throw new IllegalStateException("Byte3DBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        int numSites = values.length;
        if (numSites != myNumSites) {
            throw new IllegalArgumentException("Byte3DBuilder: setDepth: value number of sites: " + numSites + " should have: " + myNumSites);
        }
        for (int s = 0; s < myNumSites; s++) {
            set(taxon, s, values[s]);
        }
        return this;
    }

    /**
     * Add taxon and set values for all sites for that taxon.
     *
     * @param taxon taxon
     * @param values values
     *
     * @return builder
     */
    public Byte2DBuilder addTaxon(Taxon taxon, byte[] values) {
        if (myIsHDF5) {
            if (values.length != myNumSites) {
                throw new IllegalStateException("Byte3DBuilder: addTaxon: Number of sites: " + values.length + " should be: " + myNumSites);
            }
            synchronized (myHDF5Writer) {
                //HDF5Utils.writeHDF5GenotypesDepth(myHDF5Writer, taxon.getName(), values);
            }
            myNumTaxa++;
        } else {
            throw new UnsupportedOperationException();
        }
        return this;
    }

    public Byte2D build() {
        if (myIsHDF5) {
            IHDF5Reader reader = myHDF5Writer;
            myHDF5Writer = null;
            return new HDF5Byte2D(reader, mySiteScoreType);
        } else {
            SuperByteMatrix temp = myValues;
            myValues = null;
            return new MemoryByte2D(mySiteScoreType, temp);
        }
    }
}
