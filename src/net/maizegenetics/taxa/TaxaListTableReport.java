/*
 * TaxaListTableReport
 */
package net.maizegenetics.taxa;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.util.TableReport;

/**
 *
 * @author Terry Casstevens
 */
public class TaxaListTableReport implements TableReport {

    private static final String[] DEFAULT_COLUMN_HEADINGS = new String[]{"Taxa", "Name"};

    private final TaxaList myTaxaList;
    private final String[] myColumnHeadings;

    public TaxaListTableReport(TaxaList taxaList) {
        myTaxaList = taxaList;
        List<String> annotationColumns = new ArrayList<>();
        for (Taxon current : myTaxaList) {
            for (String key : current.getAnnotationAsMap().keySet()) {
                if (!annotationColumns.contains(key)) {
                    annotationColumns.add(key);
                }
            }
        }
        int totalHeadings = DEFAULT_COLUMN_HEADINGS.length + annotationColumns.size();
        myColumnHeadings = new String[totalHeadings];
        for (int i = 0; i < DEFAULT_COLUMN_HEADINGS.length; i++) {
            myColumnHeadings[i] = DEFAULT_COLUMN_HEADINGS[i];
        }
        for (int j = DEFAULT_COLUMN_HEADINGS.length; j < totalHeadings; j++) {
            myColumnHeadings[j] = annotationColumns.get(j - DEFAULT_COLUMN_HEADINGS.length);
        }
    }

    @Override
    public Object[] getTableColumnNames() {
        return myColumnHeadings;
    }

    @Override
    public String getTableTitle() {
        return "Taxa List";
    }

    @Override
    public int getColumnCount() {
        return myColumnHeadings.length;
    }

    @Override
    public int getRowCount() {
        return myTaxaList.numberOfTaxa();
    }

    @Override
    public int getElementCount() {
        return getColumnCount() * getRowCount();
    }

    @Override
    public Object[] getRow(int row) {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public Object getValueAt(int row, int col) {
        switch (col) {
            case 0:
                return myTaxaList.get(row).getName();
            case 1:
                return myTaxaList.get(row).getName();
            default:
                String[] annotations = myTaxaList.get(row).getTextAnnotation(myColumnHeadings[col]);
                if (annotations != null) {
                    return annotations[0];
                } else {
                    return null;
                }
        }
    }

}
