/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.StringTokenizer;
import java.io.FileReader;
import java.io.File;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

// import com.sleepycat.je.log.FileReader;

import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;

/**
 * This plugin takes a GBSv2 database as input, queries for the tags
 * and their taxa distribution, and creates a tab-delimited file of tag
 * and taxa-distribution.  IT can be used for verifying data in the db.
 * 
 * @author lcj34
 *
 */
public class GetTagTaxaDistFromDBPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(GetTagSequenceFromDBPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags").build();
    private PluginParameter<String> myTagFile = new PluginParameter.Builder<String>("tg", null, String.class).guiName("Input Tag File").required(true).inFile()
    		.description("Input tag file with prefered tags").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output txt file that can be imported to Excel").build();

    public GetTagTaxaDistFromDBPlugin() {
        super(null, false);
    }

    public GetTagTaxaDistFromDBPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public GetTagTaxaDistFromDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

 
    @Override
    public DataSet processData(DataSet input) {
        TagDataWriter tdw = new TagDataSQLite(inputDB());
        // Create a reader to read the TagFile
        BufferedReader tagFileReader = null;
        try {
            BufferedWriter fileWriter = null;
            StringBuilder strB = new StringBuilder();
            TaxaList taxaList = tdw.getTaxaList(); // Used for printing taxon column headers
            strB.append("Tag");
            taxaList.stream().forEach(item -> { // column names are the taxon names
                strB.append("\t");
                strB.append(item.getName());
            });
            strB.append("\n");
            // Start to read the selected tagID from the tagFile
            tagFileReader = new BufferedReader(new FileReader(new File(inputTagFile())));
            String tempString = null;
        	int srcline = 0;
        	int i = 0;
        	// Read the header of the tagFile
        	tagFileReader.readLine();
        	// Read each line to get the total tag number
        	while ((tempString = tagFileReader.readLine()) != null) {
        		StringTokenizer token=new StringTokenizer(tempString,"\t");
        		String[] array=new String[token.countTokens()];
        		i=0;
        		while(token.hasMoreTokens()){
        			array[i]=token.nextToken();
        			i++;
        		}
        		if (Integer.parseInt(array[1])>=1)
        			srcline++;
        	}
        	// System.out.println("Total tags: "+srcline);
        	tagFileReader.close();
        	// Initialize mySelTags
        	int mySelTags[] = new int[srcline];
        	// Get each Tag ID and generate mySelTags array
        	tagFileReader = new BufferedReader(new FileReader(new File(inputTagFile())));
        	// Read the header of the tagFile
        	tagFileReader.readLine();
        	while ((tempString = tagFileReader.readLine()) != null) {
        		StringTokenizer token=new StringTokenizer(tempString,"\t");
        		String[] array=new String[token.countTokens()];
        		i=0;
        		while(token.hasMoreTokens()){
        			array[i]=token.nextToken();
        			i++;
        		}
        		if (Integer.parseInt(array[1])>=1) {
        			mySelTags[srcline] = Integer.parseInt(array[1]);
        			srcline++;
        		}
        	}
        	tagFileReader.close();
            // End of generate the mySelTags array
            Set<Tag> myTags = tdw.getTags();            		
            int tagcount = 0;
            for (Tag myTag: myTags) {
                //System.out.println(myTag.toString());
            	tagcount++;
            	boolean eff = false;
	            for (int j=0;j<mySelTags.length;j++){
	            	if (tagcount == mySelTags[j]) {
	            		eff = true;
	            		break;
	            	}
	            }
	            if (eff) {
	                // get dist for each taxa
	                TaxaDistribution tagTD = tdw.getTaxaDistribution(myTag);
	                if (tagTD == null) {
	                    System.out.println("GetTagTaxaDist: got null tagTD at tagcount " + tagcount);
	                    return null;
	                }
	                int[] depths = tagTD.depths(); // gives us the depths for each taxon
	                strB.append(myTag.sequence());
	                for (int idx = 0; idx < depths.length; idx++) {
	                    // if (idx == 1 || idx == 3 || idx == 100) {
	                    strB.append("\t"); 
	                    strB.append(depths[idx]);  // add tag depth                     
	                    // }
	                }
	                strB.append("\n"); // end of line - start next tag
	            }
            }
            try {  
                fileWriter = new BufferedWriter(new FileWriter(outputFile()));
                fileWriter.write(strB.toString());
            }
            catch(IOException e) {
                myLogger.error("Caught Exception writing to outputFile " + outputFile());
                System.out.println(e);
            }
            fileWriter.close();
            ((TagDataSQLite)tdw).close();  
            myLogger.info("TagsTaxaDistToTabDelim: Finished writing TaxaDistribution \n");
        } catch (Exception exc) {
            myLogger.error("TagsTaxaDistToTabDelim: caught error " + exc);
            exc.printStackTrace();
        }
        return null;       
    }

    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
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
    public GetTagTaxaDistFromDBPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }
    /**
     * Input preferred tag file with tags and taxa distribution
     *
     * @return Input preffered tag file
     */
    public String inputTagFile() {
        return myTagFile.value();
    }
    
    /**
     * Output tab-delimited file showing tag sequences found in DB.
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. This is tab delimited file with a list
     * of tag sequences found in the database (either a single tag
     * if the user is checking for a particular tag's presence, or
     * a list of all tags found in the db (if no user tag specified).
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public GetTagTaxaDistFromDBPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

}
