/*
 *  NucleotideGenotypeTable
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.depth.AlleleDepth;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.taxa.TaxaList;

/**
 * GenotypeTable containing DNA nucleotide data.
 *
 * @author Terry Casstevens
 */
public class NucleotideGenotypeTable extends CoreGenotypeTable {

    public NucleotideGenotypeTable(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, AlleleDepth alleleDepth) {
        super(genotype, positionList, taxaList, alleleDepth);
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(genotypeAsString(taxon, i));
        }
        return builder.toString();
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
    public boolean isIndel(int site) {
        int[][] alleles = allelesSortedByFrequency(site);
        int numAlleles = Math.min(alleles[0].length, 2);
        for (int i = 0; i < numAlleles; i++) {
            if ((alleles[0][i] == NucleotideAlignmentConstants.INSERT_ALLELE) || (alleles[0][i] == NucleotideAlignmentConstants.GAP_ALLELE)) {
                return true;
            }
        }
        return false;
    }
}
