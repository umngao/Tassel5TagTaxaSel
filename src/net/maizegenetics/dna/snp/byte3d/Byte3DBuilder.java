/*
 *  Byte3DBuilder
 */
package net.maizegenetics.dna.snp.byte3d;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

import java.util.ArrayList;
import java.util.List;

/**
 * Builder to store 3 dimensional byte encoded values.
 *
 * @author Terry Casstevens
 */
public class Byte3DBuilder {

    private List<SuperByteMatrix> myValues = null;
    private final boolean myIsHDF5;
    private IHDF5Writer myHDF5Writer = null;
    private final int myMaxNumAlleles;
    private final int myNumSites;
    private int myNumTaxa = 0;

    private Byte3DBuilder(IHDF5Writer writer, int numSites, int maxNumAlleles) {
        myIsHDF5 = true;
        myHDF5Writer = writer;
        myMaxNumAlleles = maxNumAlleles;
        myNumSites = numSites;
    }

    private Byte3DBuilder(int numSites, int maxNumAlleles) {
        myIsHDF5 = false;
        myHDF5Writer = null;
        myMaxNumAlleles = maxNumAlleles;
        myNumSites = numSites;
        myValues = new ArrayList<>();
    }

    private Byte3DBuilder(int numTaxa, int numSites, int maxNumAlleles) {
        myIsHDF5 = false;
        myHDF5Writer = null;
        myMaxNumAlleles = maxNumAlleles;
        myNumSites = numSites;
        myNumTaxa = numTaxa;
        myValues = new ArrayList<>();
        for (int i = 0; i < myNumTaxa; i++) {
            myValues.add(SuperByteMatrixBuilder.getInstance(myNumSites, myMaxNumAlleles));
        }
    }

    public static Byte3DBuilder getInstance(int numTaxa, int numSites, int maxNumAlleles) {
        return new Byte3DBuilder(numTaxa, numSites, maxNumAlleles);
    }

    public static Byte3DBuilder getNucleotideInstance(int numTaxa, int numSites) {
        return new Byte3DBuilder(numTaxa, numSites, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
    }

    /**
     * Byte3DBuilder is created and values are stored in a HDF5 file.
     *
     * @param writer
     * @param numSites
     * @return
     */
    public static Byte3DBuilder getHDF5NucleotideInstance(IHDF5Writer writer, int numSites) {
        return new Byte3DBuilder(writer, numSites, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
    }

    /**
     * Set values for the all sites and alleles for a taxon simultaneously.
     * First dimension of values is number of alleles (6 for Nucleotide) and
     * second dimension is sites.
     *
     * @param taxon Index of taxon
     * @param siteOffset site offset
     * @param values array[sites][allele] of all values
     *
     * @return builder
     */
    public Byte3DBuilder setRangeForTaxon(int taxon, int siteOffset, byte[][] values) {
        if (myIsHDF5) {
            throw new IllegalStateException("Byte3DBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        int numAlleles = values.length;
        if (numAlleles != myMaxNumAlleles) {
            throw new IllegalArgumentException("Byte3DBuilder: setDepth: value number of alleles: " + numAlleles + " should have: " + myMaxNumAlleles);
        }
        for (int a = 0; a < myMaxNumAlleles; a++) {
            for (int s = 0; s < values[0].length; s++) {
                set(taxon, s + siteOffset, (byte) a, values[a][s]);
            }
        }
        return this;
    }

    /**
     * Set value for taxon, site, and allele. Value should have already been
     * translated.
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     * @param value value
     *
     * @return builder
     */
    public Byte3DBuilder set(int taxon, int site, byte allele, byte value) {
        if (myIsHDF5) {
            throw new IllegalStateException("Byte3DBuilder: set: use addTaxon for HDF5 files.");
        }
        myValues.get(taxon).set(site, allele, value);
        return this;
    }

    /**
     * Set values for the all sites and alleles for a taxon simultaneously.
     * Values should have already been translated using Byte3DUtil. First
     * dimension of values is number of alleles (6 for Nucleotide) and second
     * dimension is sites.
     *
     * @param taxon Index of taxon
     * @param values array[allele][sites] of all values
     *
     * @return builder
     */
    public Byte3DBuilder set(int taxon, byte[][] values) {
        if (myIsHDF5) {
            throw new IllegalStateException("Byte3DBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        int numAlleles = values.length;
        if (numAlleles != myMaxNumAlleles) {
            throw new IllegalArgumentException("Byte3DBuilder: setDepth: value number of alleles: " + numAlleles + " should have: " + myMaxNumAlleles);
        }
        int numSites = values[0].length;
        if (numSites != myNumSites) {
            throw new IllegalArgumentException("Byte3DBuilder: setDepth: value number of sites: " + numSites + " should have: " + myNumSites);
        }
        for (int a = 0; a < myMaxNumAlleles; a++) {
            for (int s = 0; s < myNumSites; s++) {
                set(taxon, s, (byte) a, values[a][s]);
            }
        }
        return this;
    }

    /**
     * Add taxon and set values for all sites and alleles for that taxon. First
     * dimension of values is number of alleles (6 for Nucleotide) and second
     * dimension is sites.
     *
     * @param taxon taxon
     * @param values values
     *
     * @return builder
     */
    public Byte3DBuilder addTaxon(Taxon taxon, byte[][] values) {
        if (myIsHDF5) {
            if ((values == null) || (values.length != myMaxNumAlleles)) {
                throw new IllegalStateException("Byte3DBuilder: addTaxon: Set A, C, G, T, -, + at once");
            }
            if (values[0].length != myNumSites) {
                throw new IllegalStateException("Byte3DBuilder: addTaxon: Number of sites: " + values[0].length + " should be: " + myNumSites);
            }
            synchronized (myHDF5Writer) {
                HDF5Utils.writeHDF5GenotypesDepth(myHDF5Writer, taxon.getName(), values);
            }
            myNumTaxa++;
        } else {
            myValues.add(SuperByteMatrixBuilder.getInstance(myNumSites, myMaxNumAlleles));
            set(myNumTaxa, values);
            myNumTaxa++;
        }
        return this;
    }

    public Byte3D build() {
        if (myIsHDF5) {
            IHDF5Reader reader = myHDF5Writer;
            myHDF5Writer = null;
            return new HDF5Byte3D(reader);
        } else {
            SuperByteMatrix[] temp = new SuperByteMatrix[myValues.size()];
            temp = myValues.toArray(temp);
            myValues = null;
            return new MemoryByte3D(temp, myNumTaxa, myNumSites);
        }
    }
}
