/*
 * FilterSiteNamePlugin.java
 *
 * Created on November 2, 2011
 *
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.FilterGenotypeTable;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;

import java.awt.Frame;

import java.net.URL;

import java.util.List;

import net.maizegenetics.gui.SelectFromAvailableDialog;
import net.maizegenetics.gui.SiteNamesAvailableListModel;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class FilterSiteNamePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterSiteNamePlugin.class);
    private int[] mySitesToKeep = null;
    private String[] mySiteNamesToKeep = null;
    private String[] mySiteNamesToRemove = null;

    /**
     * Creates a new instance of FilterSiteNamePlugin
     */
    public FilterSiteNamePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List inputData = input.getDataOfType(GenotypeTable.class);
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single alignment.");
                } else {
                    myLogger.error("performFunction: Please input a single alignment.");
                }
                return null;
            }

            Datum td = processDatum((Datum) inputData.get(0), isInteractive());
            if (td == null) {
                return null;
            }

            DataSet output = new DataSet(td, this);

            fireDataSetReturned(new PluginEvent(output, FilterSiteNamePlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }
    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {

        final GenotypeTable alignment = (GenotypeTable) inDatum.getData();

        if (isInteractive) {
            SiteNamesAvailableListModel listModel = new SiteNamesAvailableListModel(alignment.positions());
            SelectFromAvailableDialog dialog = new SelectFromAvailableDialog(getParentFrame(), "Site Name Filter", listModel);
            dialog.setLocationRelativeTo(getParentFrame());
            dialog.setVisible(true);
            if (dialog.isCanceled()) {
                return null;
            }
            mySitesToKeep = dialog.getDesiredIndices();
            dialog.dispose();
        }

        GenotypeTable result = null;

        if (((mySitesToKeep != null) && (mySitesToKeep.length != 0))) {
            result = FilterGenotypeTable.getInstance(alignment, mySitesToKeep);
        } else if (((mySiteNamesToKeep != null) && (mySiteNamesToKeep.length != 0))) {
            result = FilterGenotypeTable.getInstance(alignment, mySiteNamesToKeep);
        } else if (((mySiteNamesToRemove != null) && (mySiteNamesToRemove.length != 0))) {
            result = FilterGenotypeTable.getInstanceRemoveSiteNames(alignment, mySiteNamesToRemove);
        } else {
            return null;
        }

        String theName, theComment;
        theName = inDatum.getName() + "_" + result.numberOfSites() + "_Sites";
        theComment = "Subset of " + result.numberOfSites() + " from " + alignment.numberOfSites() + " Sites\n" + inDatum.getComment();
        return new Datum(theName, result, theComment);

    }

    public int[] getSitesToKeep() {
        return mySitesToKeep;
    }

    public void setSitesToKeep(int[] sitesToKeep) {
        mySitesToKeep = sitesToKeep;
        validItemsSet();
    }

    public String[] getSiteNamesToKeep() {
        return mySiteNamesToKeep;
    }

    public void setSiteNamesToKeep(String[] sitesToKeep) {
        mySiteNamesToKeep = sitesToKeep;
        validItemsSet();
    }

    public String[] getSiteNamesToRemove() {
        return mySiteNamesToRemove;
    }

    public void setSiteNamesToRemove(String[] sitesToRemove) {
        mySiteNamesToRemove = sitesToRemove;
        validItemsSet();
    }

    private void validItemsSet() {

        int count = 0;
        if ((mySitesToKeep != null) && (mySitesToKeep.length != 0)) {
            count++;
        }
        if ((mySiteNamesToKeep != null) && (mySiteNamesToKeep.length != 0)) {
            count++;
        }
        if ((mySiteNamesToRemove != null) && (mySiteNamesToRemove.length != 0)) {
            count++;
        }

        if (count > 1) {
            throw new IllegalStateException("FilterSiteNamePlugin: validItemsSet: Can only set one of the following: sites to keep, site names to keep, or site names to remove.");
        }

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterSiteNamePlugin.class.getResource("/net/maizegenetics/analysis/images/Filter.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Site Names";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Select Site Names Within Dataset";
    }
}
