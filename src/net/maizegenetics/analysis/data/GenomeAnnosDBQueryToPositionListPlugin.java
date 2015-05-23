/*
 *  GenomeAnnosDBQueryToPositionListPlugin
 * 
 *  Created on May 22, 2015
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.sql.*;
import java.util.ArrayList;
import java.util.Properties;
import java.util.stream.Collectors;
import javax.swing.*;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Jeff Glaubitz
 */
public class GenomeAnnosDBQueryToPositionListPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(GenomeAnnosDBQueryToPositionListPlugin.class);
    private static String errorMessage;

    private PluginParameter<String> connConfigFile = new PluginParameter.Builder<>("cf", null, String.class)
        .required(true)
        .inFile()
        .guiName("DB config file")
        .description("DB connection config file")
        .build();

    private PluginParameter<String> queryFile = new PluginParameter.Builder<>("qf", null, String.class)
        .required(true)
        .inFile()
        .guiName("Query file")
        .description("Query file")
        .build();

    public GenomeAnnosDBQueryToPositionListPlugin() {
        super(null, false);
    }

    public GenomeAnnosDBQueryToPositionListPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        try {
            
            myLogger.info("\nbadQueryWordsREGEX:\n"+getBadWordsRegex()+"\n\n");

            String query = readQueryFromFile();
            if (query == null) {
                complain("\nA problem occurred while reading the query file\n");
                return null;
            }

            Connection conn = connectToDB();
            if (conn == null) {
                complain("Could not connect to DB");
                return null;
            }

            PositionList posits = getPositionListBasedOnQuery(conn, query);
            if (posits == null) {
                complain("\nA problem occurred while executing the query and retreiving a PositionList\n");
                return null;
            }

            String name = Utils.getFilename(queryFile());
            Datum outputDatum = new Datum(name + "_PositionList", posits, "Position List from " + name);
            DataSet output = new DataSet(outputDatum, this);

            return output;

        } finally {
            fireProgress(100);
        }

    }
    
    private void complain(String altErrorMsg) {
        if (errorMessage == null) {
            errorMessage = altErrorMsg;
        }
        if (isInteractive()) {
            JOptionPane.showMessageDialog(getParentFrame(), errorMessage);
        } else {
            throw new IllegalStateException();
        }
    }
    
    private String readQueryFromFile() {
        BufferedReader reader = Utils.getBufferedReader(queryFile());
        String query = reader.lines().collect(Collectors.joining("\n"));
        myLogger.info("\n\nExecuting query:\n" + query + "\n\n");
        return checkQuery(query);
    }
    
    public static String checkQuery(String query) {
        String goodQuery = query.trim();
        if (!goodQuery.toUpperCase().startsWith("SELECT ")) {
            errorMessage = "\nThe supplied query must begin with \"SELECT \" (case insensitive)\n";
            myLogger.error(errorMessage);
            return null;
        }
        if (!goodQuery.toLowerCase().contains("chr")) {
            errorMessage = "\nThe supplied query must contain a output column labelled \"chr\" (lower case)\n";
            myLogger.error(errorMessage);
            return null;
        }
        if (!goodQuery.toLowerCase().contains("position")) {
            errorMessage = "\nThe supplied query must contain a output column labelled \"position\" (lower case)\n";
            myLogger.error(errorMessage);
            return null;
        }
        // ToDo: this regex-based check is not working right
        if (goodQuery.toUpperCase().matches(getBadWordsRegex())) {
            errorMessage = getBadWordsErrorMessage();
            myLogger.error(errorMessage);
            return null;
        }
        if (goodQuery.toUpperCase().contains(" DROP ")) {
            errorMessage = "\nThe supplied query must not contain the bare word \" DROP \" (case insensitive)\n";
            myLogger.error(errorMessage);
            return null;
        }
        return goodQuery;
    }
    
    private static ArrayList<String> getBadWordsArrayList() {
        ArrayList<String> badWords = new ArrayList();
        badWords.add("ALTER");
        badWords.add("COPY");
        badWords.add("CREATE");
        badWords.add("DELETE");
        badWords.add("DROP");
        badWords.add("INSERT");
        badWords.add("TRUNCATE");
        badWords.add("UPDATE");
        return badWords;
    }

    private static String getBadWordsRegex() {
        StringBuilder badWordsRegex = new StringBuilder(".*(\\s");
        badWordsRegex.append(getBadWordsArrayList().stream().collect(Collectors.joining("\\s|\\s")));
        badWordsRegex.append("\\s).*");
        return badWordsRegex.toString();
    }

    private static String getBadWordsErrorMessage() {
        StringBuilder badWordsErrMsg = new StringBuilder("Your query should not contain any of the following bare words (case insensitive):\n\t");
        badWordsErrMsg.append(getBadWordsArrayList().stream().collect(Collectors.joining("\n\t")));
        badWordsErrMsg.append("\n");
        return badWordsErrMsg.toString();
    }

    private Connection connectToDB() {
		Properties props = new Properties();
        try {
            BufferedInputStream inputStream = new BufferedInputStream(new FileInputStream(connConfigFile()));
			props.load(inputStream);
        } catch (IOException e) {
            errorMessage = "Problem reading DB connection config file (" + connConfigFile() + "):\n\t" + e;
            myLogger.error(errorMessage);
            return null;
        }
        String user = props.getProperty("user");
        String password = props.getProperty("password");
        String dbName = props.getProperty("DB");
        if (user == null) {
            errorMessage = "ERROR: Please provide a line with the user name (user=<userName>) in the DB connection config file (" + connConfigFile() + ")";
            myLogger.error(errorMessage);
            return null;
        }
        if (password == null) {
            errorMessage = "ERROR: Please provide a line with the password (password=<yourPassword>) in the DB connection config file (" + connConfigFile() + ")";
            myLogger.error(errorMessage);
            return null;
        }
        if (dbName == null) {
            errorMessage = "ERROR: Please provide a line with the DB name (DB=<dbName>) in the DB connection config file (" + connConfigFile() + ")";
            myLogger.error(errorMessage);
            return null;
        }
        return connectToDatabaseOrDie(dbName, props);
    }

    private Connection connectToDatabaseOrDie(String dbName, Properties props) {
        Connection conn = null;
        String url = "not connected yet";
        try {
            Class.forName("org.postgresql.Driver");
            url = "jdbc:postgresql://localhost/" + dbName;
            conn = DriverManager.getConnection(url, props);
        } catch (ClassNotFoundException e) {
            errorMessage = e.getMessage();
            myLogger.error(errorMessage);
            return null;
        } catch (SQLException e) {
            errorMessage = e.getMessage();
            myLogger.error(errorMessage);
            return null;
        }
        myLogger.info("\nUsing DB:  " + url + "\n");
        return conn;
    }

    private PositionList getPositionListBasedOnQuery(Connection conn, String query) {
        PositionListBuilder plb = new PositionListBuilder();
        ResultSet rs = executePostgreSQLQuery(conn, query);
        try {
            ResultSetMetaData rsmd = rs.getMetaData();
            String[] field = new String[rsmd.getColumnCount()+1];
            int[] type = new int[rsmd.getColumnCount()+1];
            for (int col = 1; col <= rsmd.getColumnCount(); col++) {
                field[col] = rsmd.getColumnLabel(col);
                type[col] = rsmd.getColumnType(col);
            }
            while (rs.next()) {
                GeneralPosition.Builder glb = new GeneralPosition.Builder(new Chromosome(rs.getString("chr")), rs.getInt("position"));
                for (int col = 1; col <= rsmd.getColumnCount(); col++) {
                    if (!field[col].equals("chr") && !field[col].equals("position")) {
                        if (type[col] == java.sql.Types.REAL || type[col] == java.sql.Types.FLOAT) {
                            glb.addAnno(field[col], rs.getFloat(col));
                        } else if (type[col] == java.sql.Types.DOUBLE) {
                            glb.addAnno(field[col], rs.getDouble(col));
                        } else {
                            glb.addAnno(field[col], rs.getString(col));
                        }
                    }
                }
                plb.add(glb.build());
            }
            rs.close();
        } catch (SQLException se) {
            errorMessage = se.getMessage();
            myLogger.error("\n"+errorMessage+"\n");
            return null;
        }
        return plb.build();
    }

    private ResultSet executePostgreSQLQuery(Connection conn, String query) {
        try {
            return conn.createStatement().executeQuery(query);
        } catch (SQLException e) {
            errorMessage = "\n\nProblem executing query (" + query + "):\n\t" + e;
            myLogger.error(errorMessage);
            return null;
        }
    }

    @Override
    public String getToolTipText() {
        return "Get a PositionList from a query of a genome annotations DB";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = GenomeAnnosDBQueryToPositionListPlugin.class.getResource("/net/maizegenetics/analysis/images/lists.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Get a PositionList from a query of a genome annotations DB";
    }
    
//    public static void main(String[] args) {
//        GeneratePluginCode.generate(GenomeAnnosDBQueryToPositionListPlugin.class);
//    }
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GenomeAnnosDBQueryToPositionListPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public PositionList runPlugin(DataSet input) {
        return (PositionList) performFunction(input).getData(0).getData();
    }

    /**
     * DB connection config file
     *
     * @return DB config file
     */
    public String connConfigFile() {
        return connConfigFile.value();
    }

    /**
     * Set DB config file. DB connection config file
     *
     * @param value DB config file
     *
     * @return this plugin
     */
    public GenomeAnnosDBQueryToPositionListPlugin connConfigFile(String value) {
        connConfigFile = new PluginParameter<>(connConfigFile, value);
        return this;
    }

    /**
     * Query file
     *
     * @return Query file
     */
    public String queryFile() {
        return queryFile.value();
    }

    /**
     * Set Query file. Query file
     *
     * @param value Query file
     *
     * @return this plugin
     */
    public GenomeAnnosDBQueryToPositionListPlugin queryFile(String value) {
        queryFile = new PluginParameter<>(queryFile, value);
        return this;
    }}
