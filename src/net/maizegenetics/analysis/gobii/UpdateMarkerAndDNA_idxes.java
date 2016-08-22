/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.io.BufferedWriter;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;

import org.apache.log4j.Logger;

import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 *
 */
public class UpdateMarkerAndDNA_idxes {
    private static final Logger myLogger = Logger.getLogger(UpdateMarkerAndDNA_idxes.class);
    public static void createIdxValues(String configFile, String outputDir, int datasetID, int platformID) {
        // connect to db
        // what should these files be called ??
        String dnarunFile = outputDir + "DS_" + datasetID + ".sh5i"; // the outputDir should include the DS_<dataset_id>
        String markerFile = outputDir + "DS_" + datasetID + ".mh5i";

        //  process the input data file
        try {
            BufferedWriter writerRunID = Utils.getBufferedWriter(dnarunFile);
            BufferedWriter writerMarkerID = Utils.getBufferedWriter(markerFile);

            long time=System.nanoTime();
            // Connect to db
            Connection dbConnection = GOBIIDbUtils.connectToDB(configFile);
            if (dbConnection == null) {
                throw new IllegalStateException("UpdateMarkerAndDNA_idxes: Problem connecting to database.");
            }
            
//            // create dnarun_id query:
//            // select dnarun_id from dataset_dnarun, dataset where dataset.name = datasetName and dataset.dataset_id = dataset_dnarun.dataset_id;
//            StringBuilder builder = new StringBuilder();
//            builder.append("select dnarun_id from dataset_dnarun, dataset where dataset.name = '");
//            builder.append(datasetName());
//            builder.append("'");
//            builder.append(" and dataset.dataset_id = dataset_dnarun.dataset_id;");
//            
//            String query = builder.toString();
//            myLogger.info("processData: query statement: " + query);
//            
//            System.out.println("UpdateMarkerAndDNA_idxes: execute query: " + query);
//            dbConnection.setAutoCommit(false); // required for Cursor processing (fetchSize)
//            Statement st = dbConnection.createStatement();
//            st.setFetchSize(100000); // should return results in batches
//            ResultSet rs = st.executeQuery(query);
//           // ResultSet rs = dbConnection.createStatement().executeQuery(query);   
//            // GOBII IFL scripts need a header in the dnarun file (but not in
//            // the .variant or .marker files !!
//            writerRunID.write("marker_name\tm_name\n"); // header
//            while (rs.next()) {
//                int dnarun_id = rs.getInt("dnarun_id");
//                writerRunID.write(String.valueOf(dnarun_id));
//                writerRunID.write("\n");
//            }
//            System.out.printf("TotalTime for dnarun_name query %g sec%n", (double) (System.nanoTime() - time) / 1e9);
            
            // create marker query:
            // select marker_id from dataset_marker, dataset where dataset.name = datasetName and dataset.dataset_id = dataset_marker.dataset_id
            StringBuilder builder = new StringBuilder();
            builder.append("select name from marker, dataset_marker where marker.marker_id=dataset_marker.marker_id and dataset_marker.dataset_id='");
            builder.append(datasetID);
            builder.append("';");
            
            String query = builder.toString();
            myLogger.info("processData: query statement: " + query);
            
            System.out.println("UpdateMarkerAndDNA_idxes: execute query: " + query);
            Statement st = dbConnection.createStatement();
            st.setFetchSize(100000); // shouldn't need to set this again
            ResultSet rs = st.executeQuery(query);
           // rs = dbConnection.createStatement().executeQuery(query); 
            writerMarkerID.write("marker_name\tm_name\tplatform_id\n"); // header
            while (rs.next()) {
                String marker_n = rs.getString("name");
                writerMarkerID.write(marker_n);
                writerMarkerID.write("\t");
                writerMarkerID.write(marker_n);
                writerMarkerID.write("\t");
                writerMarkerID.write(Integer.toString(platformID));
                writerMarkerID.write("\n");
            }
            
            writerRunID.close();
            writerMarkerID.close();
            System.out.printf("TotalTime for marker_name query: %g sec%n", (double) (System.nanoTime() - time) / 1e9);

        } catch (Exception exc) {
            System.out.println("UpdateMarkerAndDNA_idxes:  caught exception processing writing files");
            exc.printStackTrace();
        }
        System.out.println("\nFiles written to " + dnarunFile + " and " + markerFile);
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        String configFile = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/dbConfigFile_maize2.txt";
        //String datasetName = "ZeaGBSv27impV5_20160209_AGPv2_282";
        String outputDir = "/Users/lcj34/notes_files/gobiiANDBms/gobii_loading/update_idxes/";
        int datasetID = 2;
        int platformID = 2;

        createIdxValues(configFile,outputDir,datasetID,platformID); // first do the marker file
        
    }
}
