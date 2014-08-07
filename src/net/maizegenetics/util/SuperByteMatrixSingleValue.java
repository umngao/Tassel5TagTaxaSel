/*
 *  SuperByteMatrixSingleValue
 * 
 *  Created on Aug 7, 2014
 */
package net.maizegenetics.util;

import java.util.Arrays;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixSingleValue implements SuperByteMatrix {

    private byte myData;
    private final int myNumRows;
    private final int myNumColumns;

    SuperByteMatrixSingleValue(int rows, int columns, byte value) {
        myNumRows = rows;
        myNumColumns = columns;
        myData = value;
    }

    @Override
    public byte get(int row, int column) {
        return myData;
    }

    @Override
    public byte[] getAllColumns(int row) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingleValue: getAllColumns: row: " + row);
        }

        byte[] result = new byte[myNumColumns];
        Arrays.fill(result, myData);
        return result;

    }

    @Override
    public byte[] getColumnRange(int row, int start, int end) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingleValue: getColumnRange: row: " + row);
        }

        if ((start < 0) || (start >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingleValue: getColumnRange: start: " + start);
        }

        if ((end < 0) || (end >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingleValue: getColumnRange: end: " + end);
        }

        if (end < start) {
            throw new IllegalArgumentException("SuperByteMatrixSingleValue: getColumnRange: end: " + end + " less than start: " + start);
        }

        int numElements = end - start;
        byte[] result = new byte[numElements];
        Arrays.fill(result, myData);
        return result;

    }

    @Override
    public byte[] getAllRows(int column) {

        if ((column < 0) || (column >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixSingleValue: getAllRows: column: " + column);
        }

        byte[] result = new byte[myNumRows];
        for (int i = 0; i < myNumRows; i++) {
            result[i] = myData;
        }
        return result;

    }

    @Override
    public int getNumRows() {
        return myNumRows;
    }

    @Override
    public int getNumColumns() {
        return myNumColumns;
    }

    @Override
    public boolean isColumnInnerLoop() {
        return true;
    }

    @Override
    public void setHetsTo(byte value) {
        if (((myData >>> 4) & 0xf) != (myData & 0xf)) {
            myData = value;
        }
    }

    @Override
    public void set(int row, int column, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void arraycopy(int row, byte[] src, int startColumn) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setAll(byte value) {
        myData = value;
    }

    @Override
    public void reorderRows(int[] newIndices) {
        // nothing to do
    }

    @Override
    public void reorderColumns(int[] newIndices) {
        // nothing to do
    }

}
