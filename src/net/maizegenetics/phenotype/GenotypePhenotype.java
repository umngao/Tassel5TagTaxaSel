package net.maizegenetics.phenotype;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.TableReport;

/**
 * This class holds phenotypes and genotypes for a set of Taxa. 
 * @author Peter Bradbury
 * 
 */

public class GenotypePhenotype implements TableReport {
	private final GenotypeTable myGenotype;
	private final Phenotype myPhenotype;
	private final String name;
	
	/**
	 * @param theGenotype	a GenotypeTable
	 * @param thePhenotype	a Phenotype
	 * @param name	the name that will be displayed for this in the GUI
	 */
	GenotypePhenotype(GenotypeTable theGenotype, Phenotype thePhenotype, String name) {
		myGenotype = theGenotype;
		myPhenotype = thePhenotype;
		this.name = name;
	}
	
	/**
	 * @return the GenotypeTable used by this object
	 */
	public GenotypeTable genotypeTable() {
		return myGenotype;
	}
	
	/**
	 * @return	the Phenotype used by this object
	 */
	public Phenotype phenotype() {
		return myPhenotype;
	}

	/**
	 * @return	true if genotypes are discrete (nucleotides), false if the genotypes are numeric
	 * The GenotypeTable backing this object may hold both types of data. The boolean indicates which is to be used in an analysis.
	 */
	public boolean areGenotypeValuesDiscrete() {
		//TODO implement
		return true;
	}
	
	/**
	 * @return	true if taxa are replicated (more observations than taxa), false otherwise
	 */
	public boolean areTaxaReplicated() {
		//TODO implement
		return true;
	}
	
	/**
	 * @param site	the site in the GenotypeTable
	 * @return	the genotypes corresponding to every row of the phenotype table as String values
	 */
	public String[] getStringGenotype(int site) {
		//TODO implement
		return null;
	}
	
	/**
	 * @param site	the site in the GenotypeTable
	 * @return	the genotypes corresponding to every row of the phenotype table as double values
	 */
	public double[] getNumericGenotype(int site) {
		//TODO implement
		return null;
	}

	/**
	 * @param site	the site in the GenotypeTable
	 * @return	a BitSet that returns true for observations that have a missing genotype for this site, false otherwise
	 */
	public BitSet missingGenotypes(int site) {
		//TODO implement
		return null;
	}
	
	//implement TableReport methods
	@Override
	public Object[] getTableColumnNames() {
		int numberOfPhenotypeColumns = myPhenotype.getColumnCount();
		Object[] colNames = new Object[numberOfPhenotypeColumns + 1];
		System.arraycopy(myPhenotype.getTableColumnNames(), 0, colNames, 0, numberOfPhenotypeColumns);
		colNames[numberOfPhenotypeColumns] = "Genotype";
		return null;
	}

	@Override
	public String getTableTitle() {
		return name;
	}

	@Override
	public int getColumnCount() {
		return 1 + myPhenotype.getColumnCount();
	}

	@Override
	public int getRowCount() {
		return myPhenotype.getRowCount();
	}

	@Override
	public int getElementCount() {
		return getColumnCount() * getRowCount();
	}

	@Override
	public Object[] getRow(int row) {
		int numberOfPhenotypeColumns = myPhenotype.getColumnCount();
		Object[] rowData = new Object[getColumnCount()];
		System.arraycopy(myPhenotype.getRow(row), 0, rowData, 0, numberOfPhenotypeColumns);
		rowData[numberOfPhenotypeColumns] = genotypeToDisplay(row);
		return rowData;
	}

	@Override
	public Object getValueAt(int row, int col) {
        int haplotypeColumn = myPhenotype.getColumnCount();
        if (col == haplotypeColumn) return genotypeToDisplay(row);
        return myPhenotype.getValueAt(row, col);
	}
	
	private String genotypeToDisplay(int row) {
		int genotypeRow = indexOfGenotype(row);
        int siteCount = Math.min(myGenotype.numberOfSites(), 10);
        StringBuilder builder = new StringBuilder();
        builder.append(myGenotype.genotypeAsStringRange(genotypeRow, 0, siteCount));
        if (myGenotype.numberOfSites() > 10) {
            builder.append("...");
        }
        return builder.toString();
	}
	
	private int indexOfGenotype(int phenotypeRow) {
		return myGenotype.taxa().indexOf(myPhenotype.taxa().get(phenotypeRow));
	}
}
