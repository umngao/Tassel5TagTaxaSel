/**
 * 
 */
package net.maizegenetics.analysis.gobii;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;

/**
 * THis is defined in tas-1098
 * 
 * This plugin takes a haplotype breakpoint file and creates intermediate files
 * for the hap_breakpoint and breakpoint_set tables.  These three tables are currently
 * proprietary tables created by Buckler lab to be added to the GOBII postgres DB.
 * 
 * THe tables are created from the create_hapBrkptTables.sql on cbsudc01 in directory
 *   /workdir/lcj34/postgresFiles/gobii_ifl_filesToLoad/gobii_hapbreakpoints
 * 
 * THe tables to be created have these entries:
 * hapbreakpoint Table:
 *   hap_breakpoint_id  int
 *   taxa(GID)          int 
 *   positionrange      int4range (start/stop stored as an integer range)
 *   donor1 (GID)       int
 *   donor2 (GID)       int
 *   breakpoint_set_id  int  (maps to breapoint_set table)
 *   
 * breakpoint_set Table:
 *   breakpoint_set_id   int
 *   name                text
 *   method              text  (method used to created breakpoint file ,e.g FILLIN)
 *   mapset_id           int
 *   
 * projection_align Table:
 *   projection_align_id  int
 *   name                 text
 *   het_resolution       text
 *   breakpoint_set_id    int
 *   dataset_id           int  (gives us the donor file)
 * 
 * The thought is that GOBII IFL scripts will operate successfully on them.  It will
 * require the tables to be previously created in the DB, and mapping files to be
 * created and stored on the db server that can be used to process these files
 * for bulk loading.
 * 
 * THis method will take the set name and use it to populate the breakpoint table.
 * Need to create the index into each table like GOBII does, with auto increment
 * and they live in the pg_catalog.
 * 
 * Restrictions:  The breakfile format must be as defined by Ed.  The first line
 * must contain 2 tab-delimited values, the first indicating the number of donors
 * and the second indicating the number of "blocks" to process.  All lines beginning
 * with a "#" are considered comment lines.
 * 
 * There are 2 mapping files required:  one for donors, and one for taxa.  in the Ames
 * example, the donor mapping file is the same file used to curate the WGS dataset 5,
 * which is named ZeaWGS_hmp321_raw_AGPv3.  The taxa mapping file is the file Cinta
 * provides with each subset, e.g. for Ames in Cornell box.
 * 
 * The breakpoint file has donors at top:  Donors are Feb-48,PHG84, PHG83, etc.
 * Their GIDs come from the WGS mapping file Cinta gave for curating that large vcf set.
 * 1210 2711
 * #Donor Haplotypes
 * 0       Feb-48
 * 1       PHG84
 * 2       PHG83
 *
 * And taxa breakpoint blocks at the bottom:
 * "taxa" is the first column.  It gets it GID from the germplasminformation_Ames_20160810.txt 
 * file Cinta provided in Cornell Box with the Ames data: Below, 12E and 37 are taxa.
 * #Block are defined chr:startPos:endPos:donor1:donor2 (-1 means no hypothesis)
 *    12E     1:299497909:299752426:958:1070  1:299777120:300064510:79:958    1:300064521:300335541:1050:1202 1:300351838:300818858:162:7
 *    37      2:233507851:234268676:1015:1015 2:234268677:234350614:820:820   2:234350639:234411604:950:950   2:234411637:234464003:191:1
 * 
 * 
 * The third table, the projection_alignment table, will be populated with
 * another plugin.  Users will create their own projection alignment analyses
 * for the breakpoint sets of their choice.
 * @author lcj34
 *
 */
public class HapBreakpoints_IFLFilePlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(PreProcessGOBIIMappingFilePlugin.class);
    
//    private PluginParameter<String> dbConfigFile= new PluginParameter.Builder<>("dbConfigFile",null,String.class).guiName("dbConfigFile").required(true)
//            .description("DB config file containing connection information to the GOBII instance").build();
    private PluginParameter<String> breakFile= new PluginParameter.Builder<>("breakFile",null,String.class).guiName("Breakpoint File").required(true)
            .description("Full path to a single file containing breakpoint blocks for projection alignment OR to a directory containing breakpoint files,\n Files must end with *.pa.txt or *.pa.txt.gz").build();
    private PluginParameter<String> setName= new PluginParameter.Builder<>("setName",null,String.class).guiName("Breakpoint Set Name").required(true)
            .description("Name to be given to this set of breakpoints.  This name will be stored in the breakpointSet table.").build();
    private PluginParameter<String> mapset = new PluginParameter.Builder<>("mapset",null,String.class).guiName("mapset").required(true)
            .description("Name of the mapset to which these breakpoints refer.  Must match an existing name in the mapset table, e.g AGPV3").build();
    private PluginParameter<String> outputDir= new PluginParameter.Builder<>("outputDir",null,String.class).guiName("Path of output directory").required(true)
            .description("Full path name of directory to which output files will be written, must end with a /").build();
    private PluginParameter<String> method= new PluginParameter.Builder<>("method",null,String.class).guiName("Breakpoint Method").required(true)
            .description("Method used to created the breakpoints, e.g. FILLIN, beagle, etc").build();
    private PluginParameter<String> donorMapFile = new PluginParameter.Builder<>("donorMapFile",null,String.class).guiName("mappingFile").required(true)
            .description("tab-delimited File containing a TaxaColumn name to be used to map breakFile donors with a GID. \nThis may be the same file that was used for adding germplasm and marker data.").build();
    private PluginParameter<String> taxaMapFile = new PluginParameter.Builder<>("taxaMapFile",null,String.class).guiName("mappingFile").required(true)
            .description("tab-delimited File containing a TaxaColumn name to be used to map breakFile taxa with a GID. \nThis may be a different file than was used for adding germplasm and marker data.").build();
    
    List<Path> infiles;
    boolean printTaxa = false;

    public HapBreakpoints_IFLFilePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public HapBreakpoints_IFLFilePlugin() {
        super(null, false);
    }

    @Override
    protected void postProcessParameters() {
 
        File out= new File(outputDir.value());
        if (!out.getParentFile().exists()) out.getParentFile().mkdirs();
        
        // create list of directory files
        File dirList = new File(breakFile());
        if (!(dirList.exists())) {
            throw new IllegalStateException("Input file or directory not found !! " + breakFile()) ;
        }
 
        if (dirList.isDirectory()) {
            String inputFileGlob="glob:*{pa.txt,pa.txt.gz}";
            infiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(breakFile.value()).toAbsolutePath());
            if (infiles.size() < 1) {
                throw new IllegalStateException("no .txt files found in input directory !!");
            }
        } else {
            infiles=new ArrayList<>();
            Path filePath= Paths.get(breakFile()).toAbsolutePath(); 
            infiles.add(filePath);
        }
        Collections.sort(infiles); 
    } 
    //This takes the mapping file, returns a hashmap of taxaname/GID - both strings
    public HashMap<String,String> createNameGidMap(String mappingFile) {
        HashMap<String,String> taxaMapArray = new HashMap<String,String>(); 
        BufferedReader taxaMapFilebr = Utils.getBufferedReader(mappingFile, 1 << 22);
        String taxaMapLine;
        try {
            int taxaIdx=-1, gidIdx=-1;
            taxaMapLine = taxaMapFilebr.readLine();
            if (taxaMapLine != null) {
                // find header columns.  We just want taxaColumn and GID
                String [] headers = taxaMapLine.split("\\t");
                int idx = 0;
                for (String header : headers) {
                    if (header.trim().toUpperCase().equals("TAXACOLUMN")) {
                        taxaIdx = idx;
                    } else if (header.trim().toUpperCase().equals("GID")) {
                        gidIdx = idx;
                    }
                    idx++;
                }
                if (taxaIdx == -1 || gidIdx == -1) {
                    System.out.println("LCJ - createNameGIDMap returning NULL!! taxaIdx: " + taxaIdx + " gidIdx: " + gidIdx);
                    return null; // error - missing needed column
                }
            }
            // Now read file and get values
            while((taxaMapLine = taxaMapFilebr.readLine() )!= null) {
                // this is the mapping file - read and pull taxaname/gid, store in map
                String[] data = taxaMapLine.split("\\t");
                String taxaname = data[taxaIdx].trim();
                if (taxaname.contains("030-1")) {
                    System.out.println("LCJ - createNameGIdMap - found taxaname 030-1");
                }
                String gid = data[gidIdx].trim();
                if (printTaxa)
                   System.out.println("   creatNameGIDMap, adding taxa:" + taxaname + " with gid " + gid);
                taxaMapArray.put(taxaname, gid);
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        System.out.println("LCJ - size of taxaMapARray: " + taxaMapArray.entrySet().size());
        return taxaMapArray;
    }
    @Override
    public DataSet processData(DataSet input) {
        // Connect to db
        try {
//            Connection dbConnection = GOBIIDbUtils.connectToDB(dbConfigFile());
//            if (dbConnection == null) {
//                throw new IllegalStateException("PreProcessGOBIIMappingFilePlugin: Problem connecting to database.");
//            }
//            // get the dataset id;
//            StringBuilder sb = new StringBuilder();
//            sb.append("select mapset_id from mapset where name = '");
//            sb.append(mapset());
//            sb.append("';");
//            ResultSet rs = dbConnection.createStatement().executeQuery(sb.toString());
//            int mapsetId = -1;
//            while (rs.next()) {
//                mapsetId = rs.getInt("mapset_id");                
//            }
//            if (mapsetId < 0) {
//                System.out.println("Could not find mapsetId from mapset name " + mapset() + " please check name and try again !!");
//                return null;
//            }
            
            // Read the taxa mapping file, create array of names/maps
            System.out.println("LCJ - calling createNameGidMap for taxaMapFile");
            HashMap<String,String> taxaNameGid = createNameGidMap(taxaMapFile());
            if (taxaNameGid == null) {
                System.out.println("ERROR - could not create taxaname-gid mapping");
                return null;
            }
            
            //printTaxa = true;
            
            System.out.println("\nLCJ - calling createNameGIDMap for donorMapFile");
            HashMap<String,String> donorNameGid = createNameGidMap(donorMapFile());
            if (donorNameGid == null) {
                System.out.println("ERROR - could not create taxaname-gid mapping");
                return null;
            }
            System.out.println("HapBreakpoints:processData begin\n");
            String hapbrkFile = outputDir() + setName() + ".hap_breakpoint"; 
            String breakSetFile = outputDir() +  setName() + ".breakpoint_set";
            BufferedWriter writerhap = Utils.getBufferedWriter(hapbrkFile);
            BufferedWriter writerbrkset = Utils.getBufferedWriter(breakSetFile);
            
            
            // Create header lines for each intermediate file
            // Create string builders for the 3 files - we'll append data, then write
            StringBuilder hapSB = new StringBuilder();
            StringBuilder brksetSB = new StringBuilder();
            
            // write header lines: NOTE: breakpoint_set_name is for IFL to grab breakpoint_set_id value
            // this needs to appear in a hap_breakpoint.nmap file
            // also need to define the duplicate files
            hapSB.append("name\ttaxa\tpositionRange\tdonor1\tdonor2\tbreakpoint_set_name\n");
            writerhap.write(hapSB.toString());
            hapSB.setLength(0);
            
            // breakpoint_set table intermediate file. NOTE: mapset_name should be used by IFL to grab mapset_id
            // this means we don't really need to connect to DB - let IFL do it.
            // This also means don't need dbconfig parameter.  Remove it when verify is good without.
            brksetSB.append("name\tmethod\tmapset_name\n"); // header line
            brksetSB.append(setName()); // begin data line
            brksetSB.append("\t");
            brksetSB.append(method());
            brksetSB.append("\t");
            brksetSB.append(mapset());
            brksetSB.append("\n");           
            writerbrkset.write(brksetSB.toString()); // file has header and just this one line
            writerbrkset.close();
            
            BufferedReader breakerbr = null;
            String[] donorList = null;
            System.out.println("processData: number of infiles: " + infiles.size());
            for (Path filepath: infiles) {
                String curFile = filepath.toString();
                System.out.println("LCJ - processing file: " + curFile);
                breakerbr = Utils.getBufferedReader(curFile, 1 << 22);
                                                                                                                                                                                                 
                String mline = breakerbr.readLine(); // initial line containing num taxa(donor) and num blocks
                int numDonors = Integer.parseInt(mline.substring(0, mline.indexOf("\t")).trim()); // get number of donors from file header line
                int numBlocks = Integer.parseInt(mline.substring(mline.indexOf("\t")+1).trim());
                // Create list of donor GIDs. - should be same list for all files
                // Or will it be?  I think we should create the donor list from each
                // file as each file has a list
                // Process:  Read each success line, up to then number of lines 
                // the file said were there.  Grab the value after the tab (the donor nam)
                // - search for this name in the donorHashMap.  Store the donor GID
                // into the donorList array.  So the number in the "block" for the donor
                // is the index into the array, and this gives you the GID.
                System.out.println(" ... creatingDonorLists ..., number of donors is " + numDonors);
                donorList = new String[numDonors];
                boolean allFound = true;
                for (int idx = 0;idx < numDonors;idx++){
                    String donorLine = breakerbr.readLine();
                    if (donorLine.startsWith("#")) continue; // skip comment line
                    String[] donorInfo = donorLine.split("\\t");
                    String donorGid = donorNameGid.get(donorInfo[1].trim());
                    if (donorGid == null) {
                        System.out.println("LCJ - no donorGId found in donorMap for donor " + donorInfo[1]);
                        allFound = false;
                    } else {
                        donorList[idx] = donorGid;
                    }
                    donorList[idx] = donorGid;                 
                }
                if (!allFound) return null;

                System.out.println("Processing the breakpoint file ...");
                // Continue processing breakpoint file.  We processed all the donors, now process each block line
                while ( (mline=breakerbr.readLine() )!= null) {
                    hapSB.setLength(0); // reset line to 0
                    if (mline.startsWith("#")) continue;  // skip comment lines
                    // Each line contains a taxa, followed by a "block"
                    // Blocks are defined as chr:startPos:endPos:donor1:donor2 (-1 means no hypothesis)
                    // taxa and blocks are tab separated
                    String tokens[] = mline.split("\t");
                    String gidForTaxa = "-1";
                    String taxa = tokens[0].trim();
                    if (taxaNameGid.containsKey(taxa)) {
                        gidForTaxa = taxaNameGid.get(taxa);
                    }
                    if (gidForTaxa.equals("-1")) {
                        System.out.println("ERROR: gidForTaxa map does not contain taxa " + taxa + " from block list: " + mline);
                        //return null; // LCJ - want to return this error - for now I want to see if we can process anything
                        continue;
                    }
                    System.out.println("SUCCESS: found taxa " + taxa );
                    // Process blocks within each line
                    // Starting blockIdx at 1 so we pass over the taxa in 0, which was already processed
                    for (int blockIdx = 1;blockIdx < tokens.length; blockIdx++) {
                        hapSB.append(gidForTaxa); // add taxa gid
                        hapSB.append("\t");
                        // process all the blocks in this line
                        String block = tokens[blockIdx];
                        String[] blockTokens = block.split(":");
                        hapSB.append(blockTokens[0]); // add chromsome
                        hapSB.append("\t");
                        String posRange = "[" + blockTokens[1] + "," + blockTokens[2] + ")\t";
                        hapSB.append(posRange); // includes the tab
                        int donor1 = Integer.parseInt(blockTokens[3]);
                        if (donor1 >= 0 && donor1 < donorList.length) {
                            hapSB.append(donorList[donor1]);                            
                        } else {
                            System.out.println("LCJ - donor1 in block is out of range: " + block);
                            hapSB.append("-1"); // -1 for not found
                        }
                        hapSB.append("\t");
                        int donor2 = Integer.parseInt(blockTokens[4]);
                        if (donor2 >= 0 && donor2 < donorList.length) {
                            hapSB.append(donorList[donor2]);                            
                        } else {
                            System.out.println("LCJ - donor2 in block is out of range: " + block);
                            hapSB.append("-1"); // will be obvious is isn't a GID, which start at 42000000
                        }
                        hapSB.append("\t");
                        hapSB.append(setName()); // from input parameter
                        hapSB.append("\n");                                                                    
                    } 
                    // write string to file
                    writerhap.write(hapSB.toString());
                }
            }
            if (breakerbr != null) breakerbr.close(); 
            writerhap.close(); // writerbrkset closed above
        } catch (Exception exc) {
            System.out.println("processData: caught exception processing or writing files");
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

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(HapBreakpoints_IFLFilePlugin.class);
     }
     
     /**
      * Convenience method to run plugin with one return object.
      */
     // TODO: Replace <Type> with specific type.
//     public <Type> runPlugin(DataSet input) {
//         return (<Type>) performFunction(input).getData(0).getData();
//     }

//     /**
//      * DB connection config file
//      *
//      * @return dbConfigFile
//      */
//     public String dbConfigFile() {
//         return dbConfigFile.value();
//     }
//
//     /**
//      * Set dbConfigFile. DB connection config file
//      *
//      * @param value dbConfigFile
//      *
//      * @return this plugin
//      */
//     public HapBreakpoints_IFLFilePlugin dbConfigFile(String value) {
//         dbConfigFile = new PluginParameter<>(dbConfigFile, value);
//         return this;
//     }
     /**
      * Full path to file containing breakpoint blocks for
      * projection alignment
      *
      * @return Breakpoint File
      */
     public String breakFile() {
         return breakFile.value();
     }

     /**
      * Set Breakpoint File. Full path to file containing breakpoint
      * blocks for projection alignment
      *
      * @param value Breakpoint File
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin breakFile(String value) {
         breakFile = new PluginParameter<>(breakFile, value);
         return this;
     }

     /**
      * Name to be given to this set of breakpoints.  This
      * name will be stored in the breakpointSet table.
      *
      * @return Breakpoint Set Name
      */
     public String setName() {
         return setName.value();
     }

     /**
      * Set Breakpoint Set Name. Name to be given to this set
      * of breakpoints.  This name will be stored in the breakpointSet
      * table.
      *
      * @param value Breakpoint Set Name
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin setName(String value) {
         setName = new PluginParameter<>(setName, value);
         return this;
     }

     /**
      * Name of the mapset to which these breakpoints refer.
      *  Must match an existing name in the mapset table, e.g
      * AGPV3
      *
      * @return mapset
      */
     public String mapset() {
         return mapset.value();
     }

     /**
      * Set mapset. Name of the mapset to which these breakpoints
      * refer.  Must match an existing name in the mapset table,
      * e.g AGPV3
      *
      * @param value mapset
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin mapset(String value) {
         mapset = new PluginParameter<>(mapset, value);
         return this;
     }

     /**
      * Full path name of directory to which output files will
      * be written, must end with a /
      *
      * @return Path of output directory
      */
     public String outputDir() {
         return outputDir.value();
     }

     /**
      * Set Path of output directory. Full path name of directory
      * to which output files will be written, must end with
      * a /
      *
      * @param value Path of output directory
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin outputDir(String value) {
         outputDir = new PluginParameter<>(outputDir, value);
         return this;
     }

     /**
      * Method used to created the breakpoints, e.g. FILLIN,
      * beagle, etc
      *
      * @return Breakpoint Method
      */
     public String method() {
         return method.value();
     }

     /**
      * Set Breakpoint Method. Method used to created the breakpoints,
      * e.g. FILLIN, beagle, etc
      *
      * @param value Breakpoint Method
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin method(String value) {
         method = new PluginParameter<>(method, value);
         return this;
     }
     /**
      * tab-delimited File containing a TaxaColumn name to
      * be used to map breakFile donors with a GID. 
      * This may be the same file that was used for adding
      * germplasm and marker data.
      *
      * @return mappingFile
      */
     public String donorMapFile() {
         return donorMapFile.value();
     }

     /**
      * Set mappingFile. tab-delimited File containing a TaxaColumn
      * name to be used to map breakFile donors with a GID.
      * 
      * This may be the same file that was used for adding
      * germplasm and marker data.
      *
      * @param value mappingFile
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin donorMapFile(String value) {
         donorMapFile = new PluginParameter<>(donorMapFile, value);
         return this;
     }

     /**
      * tab-delimited File containing a TaxaColumn name to
      * be used to map breakFile taxa with a GID. 
      * This may be a different file than was used for adding
      * germplasm and marker data.
      *
      * @return mappingFile
      */
     public String taxaMapFile() {
         return taxaMapFile.value();
     }

     /**
      * Set mappingFile. tab-delimited File containing a TaxaColumn
      * name to be used to map breakFile taxa with a GID. 
      * This may be a different file than was used for adding
      * germplasm and marker data.
      *
      * @param value mappingFile
      *
      * @return this plugin
      */
     public HapBreakpoints_IFLFilePlugin taxaMapFile(String value) {
         taxaMapFile = new PluginParameter<>(taxaMapFile, value);
         return this;
     }


}
