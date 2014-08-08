package net.maizegenetics.trait;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.taxa.IdGroupUtils;
import net.maizegenetics.taxa.TaxaList;

/**
 * @author terry
 */
public class MarkerPhenotype implements TableReport {

    private GenotypeTable myAlignment;  //alignment to hold markers
    private Phenotype myPhenotype;  //alignment to hold characters

    /**
     * Constructor
     *
     * @param alignment alignment
     * @param phenotype phenotype
     */
    private MarkerPhenotype(GenotypeTable alignment, Phenotype phenotype) {
        myAlignment = alignment;
        myPhenotype = phenotype;
    }

    public static MarkerPhenotype getInstance(GenotypeTable aa, Phenotype ca, boolean union) {
        TaxaList idGroup = getIdGroup(aa.taxa(), ca.getTaxa(), union);
        GenotypeTable align = FilterGenotypeTable.getInstance(aa, idGroup);
        Phenotype phenotype = FilterPhenotype.getInstance(ca, idGroup, null);
        return new MarkerPhenotype(align, phenotype);
    }

    public static MarkerPhenotype getInstance(MarkerPhenotype aac, TaxaList group) {
        GenotypeTable aa = FilterGenotypeTable.getInstance(aac.getAlignment(), group);
        Phenotype ca = FilterPhenotype.getInstance(aac.getPhenotype(), group, null);
        return new MarkerPhenotype(aa, ca);
    }

    private static TaxaList getIdGroup(TaxaList group1, TaxaList group2, boolean union) {
        if (union) {
            return IdGroupUtils.getAllIds(group1, group2);
        } else {
            return IdGroupUtils.getCommonIds(group1, group2);
        }
    }

    public String toString() {

        StringBuilder builder = new StringBuilder();

        builder.append("  ");
        builder.append(myAlignment.numberOfTaxa());
        builder.append(" ");
        builder.append(myPhenotype.getNumberOfTraits());
        builder.append("Taxa/Trait\t");
        for (int j = 0; j < myPhenotype.getNumberOfTraits(); j++) {
            builder.append(myPhenotype.getTrait(j).getName());
            builder.append("\t");
        }
        builder.append("\n");
        for (int i = 0; i < myAlignment.numberOfTaxa(); i++) {
            builder.append(myAlignment.taxaName(i));
            builder.append("\t");
            for (int j = 0; j < myPhenotype.getNumberOfTraits(); j++) {
                builder.append(myPhenotype.getData(i, j));
                builder.append("\t");
            }
            builder.append(myAlignment.genotypeAsStringRow(i));
        }
        return builder.toString();
    }

    public GenotypeTable getAlignment() {
        return myAlignment;
    }

    public Phenotype getPhenotype() {
        return myPhenotype;
    }

    public Object[] getTableColumnNames() {
        String[] basicLabels = new String[myPhenotype.getNumberOfTraits() + 2];
        basicLabels[0] = "Taxa";
        for (int c = 0; c < myPhenotype.getNumberOfTraits(); c++) {
            basicLabels[c + 1] = myPhenotype.getTrait(c).getName();
        }
        basicLabels[myPhenotype.getNumberOfTraits() + 1] = "Haplotype";
        return basicLabels;
    }

    public String getTableTitle() {
        return "Phenotypes and Genotypes";
    }

    public int getColumnCount() {
        return myPhenotype.getNumberOfTraits() + 2;
    }

    public int getRowCount() {
        return myAlignment.numberOfTaxa();
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public Object[] getRow(int row) {

        Object[] data;
        data = new String[myPhenotype.getNumberOfTraits() + 2];
        data[0] = myAlignment.taxaName(row);
        for (int c = 0; c < myPhenotype.getNumberOfTraits(); c++) {
            data[c + 1] = "" + myPhenotype.getData(row, c);
        }
        int siteCount = Math.min(myAlignment.numberOfSites(), 10);
        StringBuilder builder = new StringBuilder();
        builder.append(myAlignment.genotypeAsStringRange(row, 0, siteCount));
        if (myAlignment.numberOfSites() > 10) {
            builder.append("...");
        }
        data[myPhenotype.getNumberOfTraits() + 1] = builder.toString();
        return data;

    }

    public Object getValueAt(int row, int col) {
        int haplotypeColumn = myPhenotype.getColumnCount();
        if (col == haplotypeColumn) {
            int siteCount = Math.min(myAlignment.numberOfSites(), 10);
            StringBuilder builder = new StringBuilder();
            builder.append(myAlignment.genotypeAsStringRange(row, 0, siteCount));
            if (myAlignment.numberOfSites() > 10) {
                builder.append("...");
            }
            return builder.toString();
        }
        return myPhenotype.getValueAt(row, col);
    }
}
