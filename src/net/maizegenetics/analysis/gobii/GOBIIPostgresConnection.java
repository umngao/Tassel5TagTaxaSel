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
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
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

    public static TaxaList taxaList(Connection connection, String datasetName) {

        if (connection == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: taxaList: Must specify database connection.");
        }

        // select distinct(germplasm.name) from dataset, dataset_dnarun, dnarun,
        // dnasample, germplasm
        // where dataset.name='maize282_raw_AGPv2'
        // and dataset_dnarun.dataset_id = dataset.dataset_id
        // and dnarun.dnarun_id = dataset_dnarun.dnarun_id
        // and dnarun.dnasample_id = dnasample.dnasample_id
        // and dnasample.germplasm_id = germplasm.germplasm_id;
        StringBuilder builder = new StringBuilder();
        builder.append("select distinct(germplasm.name) from dataset, dataset_dnarun, dnarun, dnasample, germplasm ");
        builder.append("where dataset.name='");
        builder.append(datasetName);
        builder.append("'");
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

    public static PositionList positionList(Connection connection, String datasetName) {

        if (connection == null) {
            throw new IllegalArgumentException("GOBIIPostgresConnection: positionList: Must specify database connection.");
        }

        // select marker.name, marker_linkage_group.start, marker_linkage_group.stop, linkage_group.name
        // from dataset, dataset_marker, marker, marker_linkage_group, linkage_group
        // where dataset.name='maize282_raw_AGPv2'
        // and dataset.dataset_id=dataset_marker.dataset_id
        // and dataset_marker.marker_id=marker.marker_id
        // and marker.marker_id=marker_linkage_group.marker_id
        // and marker_linkage_group.linkage_group_id=linkage_group.linkage_group_id;
        StringBuilder builder = new StringBuilder();
        builder.append("select marker.name, marker_linkage_group.start, marker_linkage_group.stop, linkage_group.name ");
        builder.append("from dataset, dataset_marker, marker, marker_linkage_group, linkage_group ");
        builder.append("where dataset.name='");
        builder.append(datasetName);
        builder.append("'");
        builder.append(" and dataset.dataset_id=dataset_marker.dataset_id");
        builder.append(" and dataset_marker.marker_id=marker.marker_id");
        builder.append(" and marker.marker_id=marker_linkage_group.marker_id");
        builder.append(" and marker_linkage_group.linkage_group_id=linkage_group.linkage_group_id;");

        String query = builder.toString();
        myLogger.info("taxaList: query statement: " + query);

        try (ResultSet rs = connection.createStatement().executeQuery(query)) {
            PositionListBuilder positions = new PositionListBuilder();
            while (rs.next()) {
                int start = rs.getInt("start");
                int end = rs.getInt("stop");
                if (start != end) {
                    throw new IllegalArgumentException("GOBIIPostgresConnection: positionList: start position: " + start + " and end position: " + end + " should be the same.");
                }
                String snpName = rs.getString(1);
                // Todo - Insert Chromosome
                GeneralPosition.Builder current = new GeneralPosition.Builder(Chromosome.UNKNOWN, start);
                current.snpName(snpName);
                positions.add(current.build());
            }
            return positions.build();
        } catch (Exception se) {
            myLogger.debug(se.getMessage(), se);
            throw new IllegalStateException("GOBIIPostgresConnection: taxaList: Problem querying the database: " + se.getMessage());
        }

    }

    public static void main(String[] args) {

        String datasetName = "maize282_raw_AGPv2";

        Connection conneciton = connection("/home/tmc46/gobii_maizeifltest_props.txt");

        TaxaList taxa = taxaList(conneciton, datasetName);
        System.out.println("number of taxa: " + taxa.numberOfTaxa());

        PositionList positions = positionList(conneciton, datasetName);
        System.out.println("number of positions: " + positions.numberOfSites());

    }

}
