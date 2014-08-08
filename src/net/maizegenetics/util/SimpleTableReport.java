package net.maizegenetics.util;

import java.io.Serializable;

/**
 * Created by IntelliJ IDEA. User: ed Date: Sep 28, 2006 Time: 9:37:46 PM
 */
public class SimpleTableReport extends AbstractTableReport implements Serializable, TableReport {

    Object[][] theData;
    Object[] theColumnNames;
    Object[] theRowNames = null;
    String theName;

    //Implementation of TableReport Interface
    public SimpleTableReport(String theName, Object[] columnNames, Object[][] theData) {
        this.theData = theData;
        this.theColumnNames = columnNames;
        this.theName = theName;
    }

    public SimpleTableReport(TableReport tr) {
        theData = new Object[tr.getRowCount()][tr.getColumnCount()];
        int numRows = tr.getRowCount();
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(tr.getRow(i), 0, theData[i], 0, numRows);
        }
        theColumnNames = tr.getTableColumnNames();
        theName = tr.getTableTitle();
    }

    /**
     * Return column names for the table
     */
    public Object[] getTableColumnNames() {
        return theColumnNames;
    }

    /**
     * Return data for the table
     */
    public Object[][] getTableData() {
        return theData;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {
        return theData[row];
    }

    /**
     * Return the name for the title of the ANOVA
     */
    public String getTableTitle() {
        return theName;
    }

    public int getRowCount() {
        return theData.length;
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public int getColumnCount() {
        return theColumnNames.length;
    }

    public void setRowNames(Object[] rowNames) {
        theRowNames = rowNames;
    }

    public Object getValueAt(int row, int col) {
        return theData[row][col];
    }
}
