// TableReport.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.util;

/**
 * Interface for classes with data that can be presented in tables
 *
 * @author Ed Buckler
 */
public interface TableReport {

    /**
     * get the names of the columns
     *
     * @return columns names
     */
    public Object[] getTableColumnNames();

    /**
     * get the title of the table
     *
     * @return a String title
     */
    public String getTableTitle();

    /**
     * get the number of the columns
     *
     * @return number of columns
     */
    public int getColumnCount();

    /**
     * get the number of rows
     *
     * @return number of rows
     */
    public int getRowCount();

    /**
     * Get the total number of elements in the dataset. Elements=rowCount *
     * columnCount;
     *
     * @return columns names
     */
    public int getElementCount();

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row);

    /**
     * Returns value at given row and column.
     *
     * @param row row number
     * @param col column number
     * @return data
     */
    public Object getValueAt(int row, int col);

}
