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
 * Once I have the datasets fixed, this functionality should not be needed.
 * 
 * What it does:  Initially the marker_idx and dnarun_idx coluns of the dataset_marker
 * and dataset_dnarun tables respectively were not populated.  They are now needed and
 * are populated.  Kevin Palis created a couple scripts to handle populating these fields
 * in tables when they were missing.  These scripts live with the gobii_ifl_scripts
 * on CBSU, and are called update_marker_idx.py and update_dnarun_idx.py.  For GOBII,
 * then are in the rpository at the same level as the gobii_ifl.py scrips.
 * 
 * The file below reates an intermediate file that will be worked on by the preprocess_ifile.py
 * script.  You can also run the gobii_ifl.py script instead if you uncomment the "return" statement
 * that occurs after the preprocess_ifile.py script has been called.
 * 
 * 
 * @author lcj34
 *
 */
public class UpdateMarkerAndDNA_idxes {
    private static final Logger myLogger = Logger.getLogger(UpdateMarkerAndDNA_idxes.class);
    public static void createIdxValues(String configFile, String outputDir, int datasetID, int platformID, int experimentID) {
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
            
            // create dnarun_idx query:  here we need the dnarun.name field
//            StringBuilder builder = new StringBuilder();
//            
//            builder.append("select name from dnarun,dataset_dnarun where dataset_dnarun.dnarun_id=dnarun.dnarun_id and dataset_dnarun.dataset_id='");
//            builder.append(datasetID);
//            
//            String query = builder.toString();
//            myLogger.info("processData: query statement for dnarun: " + query);
//            
//            System.out.println("UpdateMarkerAndDNA_idxes: execute query: " + query);
//            dbConnection.setAutoCommit(false); // required for Cursor processing (fetchSize)
//            Statement st = dbConnection.createStatement();
//            st.setFetchSize(100000); // should return results in batches
//            ResultSet rs = st.executeQuery(query);
//           // ResultSet rs = dbConnection.createStatement().executeQuery(query);   
//            // GOBII IFL scripts need a header in the dnarun file (but not in
//            // the .variant or .marker files !!
//            writerRunID.write("dnarun_name\td_name\texperiment_id\n"); // header
//            while (rs.next()) {
//                String dnarun_n = rs.getString("name");
//                writerRunID.write(dnarun_n);
//                writerRunID.write("\t");
//                writerRunID.write(dnarun_n); // write name out twice as first one is converted.
//                writerRunID.write("\t");
//                writerRunID.write(Integer.toString(experimentID));
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
                writerMarkerID.write(marker_n); // write same name twice
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
        int datasetID = 4;
        int platformID = 3; // needed for marker file
        int experimentID = 1; // needed for dnarun file

        createIdxValues(configFile,outputDir,datasetID,platformID, experimentID); // first do the marker file
        
    }
}
