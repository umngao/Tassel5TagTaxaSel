/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Ordering;

/**
 * This plugin takes as input:
 *    csv file with columns chr, pos, qualityScore 
 *    dbFile:  a GBSv2 database with snppositions recorded
 *  A map of the chromosome, position and qualityScore is sent to the database
 *  where the snpposition table is updated with a qualityScore value for the
 *  specified chromosome and position.
 *  
 * @author lcj34
 *
 */
public class UpdateSNPPositionQualityPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(UpdateSNPPositionQualityPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with SNP positions stored").build();
    private PluginParameter<String> myQSFile = new PluginParameter.Builder<String>("qsFile", null, String.class).guiName("Quality Score File").required(true).inFile()
            .description("CSV file containing headers chr(String), pos(Integer) and qualityScore(Float) for filtering SNP positions from database").build();

    public UpdateSNPPositionQualityPlugin() {
        super(null, false);
    }

    public UpdateSNPPositionQualityPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public UpdateSNPPositionQualityPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        try {
            Path qsPath= Paths.get(qsFile()).toAbsolutePath();
            ListMultimap<String, Tuple<Integer, Float>> qsMap = readQualityScoreFile(qsPath.toString());
            if (qsMap == null) {
            	String errMsg = "Error reading csv input file " + qsPath.toString();
            	myLogger.error(errMsg);
            	return null;
            }
            TagDataWriter tdw=new TagDataSQLite(inputDB());
            // Write the qualityScore positions to the snpposition table
            tdw.putSNPPositionQS(qsMap);           
            ((TagDataSQLite)tdw).close();  
            myLogger.info("Finished writing quality scores file to snpposition table.\n");
        } catch (Exception exc) {
            myLogger.error("Caught error adding quality scores to the database " + exc);
            exc.printStackTrace();
        }
        return null;
    }

    /**
     * Parses a csv qualityScore file storing the chr, position and quality score values into a multimap.
     * The chr (chromosome) is the key, which may have multiple associated position/quality scores.
     *
     * @param filename: full name of csv file for processing
     * @return A multi-map with String (Chromosome) as key, and a tuple containing position and qualityScore
     */
    public static ListMultimap<String, Tuple<Integer, Float>> readQualityScoreFile(String fileName) {
        if (fileName == null) {
            return null;
        }
        ImmutableListMultimap.Builder<String, Tuple<Integer, Float>> qsMap = new ImmutableListMultimap.Builder<String, Tuple<Integer, Float>>()
                .orderKeysBy(Ordering.natural()).orderKeysBy(Ordering.natural());
        try {
            BufferedReader br = Utils.getBufferedReader(fileName, 1000000);
            String line;
            while((line=br.readLine())!=null) {
                if(line.contains("qualityScore")) continue; // skip header
                String[] myString = line.split(",");
                if (myString.length < 3) {
                	String msg = "ERROR: wrong number of fields in Quality Score File.\n" +
                			"Expecting chr,pos,qualityScore but found:\n" + line + "\n";
                	myLogger.error(msg);
                	return null;
                }
                String myChr = myString[0];
                Integer myPos = Integer.parseInt(myString[1]);
                Float myQS = Float.parseFloat(myString[2]);
                Tuple<Integer,Float> myTuple = new Tuple<Integer,Float>(myPos,myQS);
                qsMap.put(myChr, myTuple);
            }
            br.close();           
        } catch (Exception e) {
            System.err.println("Error in Reading QualityScore  File:" + fileName);
            e.printStackTrace();
        }
        return qsMap.build();
    }

    @Override
    public String getToolTipText() {
        return "Update SNP Position quality score";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Update SNP Position Quality";
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
    }

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public UpdateSNPPositionQualityPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Input quality score file to use for setting the qualityScore
     * in the database's snpposition table.
     *
     * @return QualityScore File
     */
    public String qsFile() {
        return myQSFile.value();
    }

    /**
     * Set qsFile. Input quality score file used for setting
     * the qualityScore field in the snpposition table.
     *
     * @param value Quality Score File
     *
     * @return this plugin
     */
    public UpdateSNPPositionQualityPlugin qsFile(String value) {
        myQSFile = new PluginParameter<>(myQSFile, value);
        return this;
    }
}
