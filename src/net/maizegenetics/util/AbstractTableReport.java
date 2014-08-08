/*
 * AbstractTableReport
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractTableReport implements TableReport {

    private int currentRowNumber = -1;
    private Object[] currentRow = null;

    public Object getValueAt(int row, int col) {
        if (row != currentRowNumber) {
            currentRowNumber = row;
            currentRow = getRow(row);
        }
        return currentRow[col];
    }

}
