/*
 *  GOBIIPlugin
 * 
 *  Created on July 25, 2016
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.sql.Connection;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GOBIIGenotypeCallTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.TaxaList;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class GOBIIPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(GOBIIPlugin.class);

    private PluginParameter<String> myDataset = new PluginParameter.Builder<>("dataset", null, String.class)
            .required(true)
            .description("Dataset Name")
            .guiName("Dataset Name")
            .build();

    private final PluginParameter<String> myDatabaseLabel = PluginParameter.getLabelInstance("Database Properties");

    private PluginParameter<String> myDBName = new PluginParameter.Builder<>("db", null, String.class)
            .required(true)
            .guiName("Database Name")
            .description("Database Name")
            .build();

    private PluginParameter<String> myUser = new PluginParameter.Builder<>("user", null, String.class)
            .required(true)
            .description("User Name")
            .build();

    private PluginParameter<String> myPassword = new PluginParameter.Builder<>("password", "", String.class)
            .description("Password")
            .password()
            .build();

    public GOBIIPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if ((dBName() == null) || (dBName().length() == 0)) {
            dBName(TasselPrefs.getGOBIIDB());
        }
        if ((user() == null) || (user().length() == 0)) {
            user(TasselPrefs.getGOBIIUser());
        }
    }

    @Override
    protected void postProcessParameters() {
        TasselPrefs.putGOBIIDB(dBName());
        TasselPrefs.putGOBIIUser(user());
    }

    @Override
    public DataSet processData(DataSet input) {

        Connection conneciton = GOBIIPostgresConnection.connection("localhost", user(), myPassword.value(), dBName());

        String hdf5Filename = null;
        try {
            hdf5Filename = GOBIIPostgresConnection.hdf5Filename(conneciton, dataset());
        } catch (Exception e) {
            GOBIIPostgresConnection.printAvailableDatasets(conneciton);
            throw e;
        }

        TaxaList taxa = GOBIIPostgresConnection.taxaList(conneciton, dataset());
        myLogger.info("Number of Taxa: " + taxa.numberOfTaxa());

        PositionList positions = GOBIIPostgresConnection.positionList(conneciton, dataset());
        myLogger.info("Number of Positions: " + positions.numberOfSites());

        GOBIIGenotypeCallTable genotypes = GOBIIGenotypeCallTable.getInstance(taxa.numberOfTaxa(), positions.numberOfSites(), false, hdf5Filename);

        GenotypeTable genotypeTable = GenotypeTableBuilder.getInstance(genotypes, positions, taxa);
        return new DataSet(new Datum("gobii:" + dataset(), genotypeTable, null), this);

    }

    /**
     * Dataset Name
     *
     * @return Dataset Name
     */
    public String dataset() {
        return myDataset.value();
    }

    /**
     * Set Dataset Name. Dataset Name
     *
     * @param value Dataset Name
     *
     * @return this plugin
     */
    public GOBIIPlugin dataset(String value) {
        myDataset = new PluginParameter<>(myDataset, value);
        return this;
    }

    /**
     * Database Name
     *
     * @return Database Name
     */
    public String dBName() {
        return myDBName.value();
    }

    /**
     * Set Database Name. Database Name
     *
     * @param value Database Name
     *
     * @return this plugin
     */
    public GOBIIPlugin dBName(String value) {
        myDBName = new PluginParameter<>(myDBName, value);
        return this;
    }

    /**
     * User Name
     *
     * @return User
     */
    public String user() {
        return myUser.value();
    }

    /**
     * Set User. User Name
     *
     * @param value User
     *
     * @return this plugin
     */
    public GOBIIPlugin user(String value) {
        myUser = new PluginParameter<>(myUser, value);
        return this;
    }

    /**
     * Set Password. Password
     *
     * @param value Password
     *
     * @return this plugin
     */
    public GOBIIPlugin password(String value) {
        myPassword = new PluginParameter<>(myPassword, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "GOBII Plugin";
    }

    @Override
    public String getToolTipText() {
        return "GOBII Plugin";
    }

}
