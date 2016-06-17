/*
 *  GOBIIPostgresConnection
 * 
 *  Created on Jun 15, 2016
 */
package net.maizegenetics.analysis.gobii;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Properties;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class GOBIIPostgresConnection {

    private static final Logger myLogger = Logger.getLogger(GOBIIPostgresConnection.class);

    private GOBIIPostgresConnection() {
        //utility
    }

    public static Connection connection(String propertiesFile) {

        Properties properties = new Properties();
        try {
            properties.load(Utils.getBufferedReader(propertiesFile));
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("GOBIIPostgresConnection: connection: problem reading properties file: " + propertiesFile);
        }

        String host = properties.getProperty("host");
        String user = properties.getProperty("user");
        String password = properties.getProperty("password");
        String dbName = properties.getProperty("DB");

        Connection connection = null;
        String url = "jdbc:postgresql://localhost/" + dbName;
        try {
            Class.forName("org.postgresql.Driver");
            connection = DriverManager.getConnection(url, user, password);
        } catch (ClassNotFoundException e) {
            myLogger.error(e.getMessage());
            return null;
        } catch (SQLException e) {
            myLogger.error(e.getMessage());
            return null;
        }
        myLogger.info("Connected to database:  " + url + "\n");
        return connection;

    }

    public static TaxaList taxaList(Connection connection, int experimentID, int analysisID) {

        if (connection == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: taxaList: Must specify database connection.");
        }

        // select distinct(germplasm.name) from dataset, dataset_dnarun, dnarun,
        // dnasample, germplasm
        // where dataset.experiment_id=5 and dataset.callinganalysis_id=4
        // and dataset_dnarun.dataset_id = dataset.dataset_id
        // and dnarun.dnarun_id = dataset_dnarun.dnarun_id
        // and dnarun.dnasample_id = dnasample.dnasample_id
        // and dnasample.germplasm_id = germplasm.germplasm_id;
        StringBuilder builder = new StringBuilder();
        builder.append("select distinct(germplasm.name) from dataset, dataset_dnarun, dnarun, dnasample, germplasm ");
        builder.append("where dataset.experiment_id=");
        builder.append(experimentID);
        builder.append(" and dataset.callinganalysis_id=");
        builder.append(analysisID);
        builder.append(" and dataset_dnarun.dataset_id = dataset.dataset_id");
        builder.append(" and dnarun.dnarun_id = dataset_dnarun.dnarun_id");
        builder.append(" and dnarun.dnasample_id = dnasample.dnasample_id");
        builder.append(" and dnasample.germplasm_id = germplasm.germplasm_id;");

        String query = builder.toString();
        myLogger.info("taxaList: query statement: " + query);

        try (ResultSet rs = connection.createStatement().executeQuery(query)) {
            TaxaListBuilder taxa = new TaxaListBuilder();
            while (rs.next()) {
                Taxon current = new Taxon(rs.getString("name"));
                taxa.add(current);
            }
            return taxa.build();
        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("GOBIIPostgresConnection: taxaList: Problem querying the database: " + se.getMessage());
        }

    }

    public static void main(String[] args) {
        int experimentID = 1;
        int analysisID = 1;
        Connection conneciton = connection("/home/tmc46/gobii_maizeifltest_props.txt");
        TaxaList taxa = taxaList(conneciton, experimentID, analysisID);
        System.out.println("number of taxa: " + taxa.numberOfTaxa());
    }

}
