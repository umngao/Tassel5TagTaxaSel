package net.maizegenetics.phenotype;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.TableReport;

public class GenotypePhenotype implements TableReport {
	private final GenotypeTable myGenotype;
	private final Phenotype myPhenotype;
	private final String name;
	
	GenotypePhenotype(GenotypeTable theGenotype, Phenotype thePhenotype, String name) {
		myGenotype = theGenotype;
		myPhenotype = thePhenotype;
		this.name = name;
	}
	
	public GenotypeTable genotypeTable() {
		return myGenotype;
	}
	
	public Phenotype phenotype() {
		return myPhenotype;
	}

	public boolean areGenotypeValuesDiscrete() {
		//TODO implement
		return true;
	}
	
	public String[] getStringGenotype() {
		//TODO implement
		return null;
	}
	
	public String[] getStringGenotype(boolean[] notMissing) {
		//TODO implement
		return null;
	}
	
	public double[] getNumericGenotype() {
		//TODO implement
		return null;
	}
	
	public double[] getNumericGenotype(boolean[] notMissing) {
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
	public Object[][] getTableData() {
		getTableData(0, getRowCount() - 1);
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
	public Object[][] getTableData(int start, int end) {
		int nrows = end - start + 1;
		Object[][] tableData = new Object[nrows][];
		for (int i = 0; i < nrows; i++) tableData[i] = getRow(i + start);
		return null;
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
