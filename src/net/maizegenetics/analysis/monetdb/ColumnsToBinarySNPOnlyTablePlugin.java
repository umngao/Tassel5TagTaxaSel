/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;

/**
 * @author lcj34
 *
 */
public class ColumnsToBinarySNPOnlyTablePlugin extends AbstractPlugin {
    private PluginParameter<String> inputFile= new PluginParameter.Builder<>("inputFile",null,String.class).guiName("Input File").required(true)
            .description("Input File or Directory containing Tab-delimited files with data to add to the database. Files must be named chr01.txt, chr02.txt etc!").build();
    private PluginParameter<String> refDir= new PluginParameter.Builder<>("refDir",null,String.class).guiName("refDir").required(true)
            .description("Directory holding files that contain all the SNP physical positions, separated by chr").build();
    private PluginParameter<String> outBase= new PluginParameter.Builder<>("outBase",null,String.class).guiName("outBase").required(true)
            .description("Output directory and base filename to hold the binary files. Will make directory if neccesary").build();
    private PluginParameter<String> colsFloat= new PluginParameter.Builder<>("colNamesFloat",null,String.class).guiName("Columns keep as Real (Float)").required(false)
            .description("Comma separated list of column names to generate real binaries for").build();
    private PluginParameter<String> colsInt= new PluginParameter.Builder<>("colNamesInt",null,String.class).guiName("Columns keep as Int").required(false)
            .description("Comma separated list of column names to generate int  binaries for").build();
    private PluginParameter<String> colsShort= new PluginParameter.Builder<>("colNamesShort",null,String.class).guiName("Columns keep as Short").required(false)
            .description("Comma separated list of column names to generate short binaries for").build();
    private PluginParameter<String> colsLong= new PluginParameter.Builder<>("colNamesLong",null,String.class).guiName("Columns keep as Long").required(false)
            .description("Comma separated list of column names to generate long binaries for").build();
    private PluginParameter<String> colsChar= new PluginParameter.Builder<>("colNamesChar",null,String.class).guiName("Columns keep as Char String").required(false)
            .description("Comma separated list of column names to generate character/text binaries for").build();
    private PluginParameter<String> colsByte= new PluginParameter.Builder<>("colNamesByte",null,String.class).guiName("Columns keep as Byte").required(false)
            .description("Comma separated list of column names to generate byte binaries for").build();
    private PluginParameter<String> colsLog10= new PluginParameter.Builder<>("colNamesLog10",null,String.class).guiName("Columns to Keep and transform -log10").required(false)
            .description("Comma separated list of column names to first transform using -log10 then generate binaries for").build();
    private PluginParameter<Boolean> range= new PluginParameter.Builder<>("range",false,Boolean.class).guiName("Range information?").required(false)
            .description("Columns for range data. If true, will look for 'start' and 'end' (inclusive exclusive) or 'first' 'last' (inclusive inclusive) instead of 'Pos' ").build();
    private PluginParameter<Boolean> negToZero= new PluginParameter.Builder<>("negToZero",false,Boolean.class).guiName("Negative floats to zero?").required(false)
            .description("Will set negative float denoted column values to zero (otherwise they are set to Float.MIN").build();
    private PluginParameter<Boolean> missToZero= new PluginParameter.Builder<>("missToZero",false,Boolean.class).guiName("Missing float values to zero?").required(false)
            .description("Will set missing float denoted column values to zero (otherwise they are set to Float.MIN").build();
    private PluginParameter<Boolean> oneBased= new PluginParameter.Builder<>("oneBased",true,Boolean.class).guiName("positions are 1 based?").required(false)
            .description("Will assume all positions and ranges are 1-based unless this value is set to false").build();
    
    private HashMap<Integer,String> log10HM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> log10Writers= null;
    private HashMap<Integer,String> realHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> realWriters= null;
    private HashMap<Integer,String> shortHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> shortWriters= null;
    private HashMap<Integer,String> longHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> longWriters= null;
    ///PUT CHARACTER BACK (figure out what monetDB wants)
    private HashMap<Integer,String> charHM= null; //Strings written as Strings, because that's what MonetDB likes
    private HashMap<Integer,DataOutputStream> charWriters= null;
    private HashMap<Integer,String> intHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> intWriters= null;
    private HashMap<Integer,String> byteHM= null;
    private HashMap<Integer,LittleEndianDataOutputStream> byteWriters= null;
    private int chrCol= -1;
    private String chrName= "CHR";
    private int posCol= -1;
    private int totalLines;
    private int totalZeroLines;
    private int totalNonZeroLines;
    private int linesForChr;
    private FileInputStream fis= null;
    private Scanner scanner= null;
    private boolean inclusive= false;
    private int startCol= -1;
    private int endCol= -1;
    
    List<Path> infiles;
//    private boolean firstChar= true;
    
    public ColumnsToBinarySNPOnlyTablePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    public ColumnsToBinarySNPOnlyTablePlugin() {
        super(null, false);
    }
    
    @Override
    protected void preProcessParameters(DataSet input) {
    }
    
    @Override
    protected void postProcessParameters() {
 
        File out= new File(outBase.value());
        if (!out.getParentFile().exists()) out.getParentFile().mkdirs();
        
        // create list of directory files
        File dirList = new File(inputFile());
        if (!(dirList.exists())) {
            throw new IllegalStateException("Input file or directory not found !!") ;
        }
 
        if (dirList.isDirectory()) {
            String inputFileGlob="glob:*{.txt}";
            infiles= DirectoryCrawler.listPaths(inputFileGlob, Paths.get(inputFile.value()).toAbsolutePath());
            if (infiles.size() < 1) {
                throw new IllegalStateException("no .txt files found in input directory !!");
            }
        } else {
            infiles=new ArrayList<>();
            Path filePath= Paths.get(inputFile()).toAbsolutePath(); 
            infiles.add(filePath);
        }
        Collections.sort(infiles); // files should be chr01.txt ...chr10.txt so they will sort correctly
        // all chroms should have same file headers - get columns from first file
        // This is a requirement to have headers in the files.
        FindColumns(infiles.get(0).toString());
    }

    @Override
    public DataSet processData(DataSet input) {
        long time= System.nanoTime();
        Initialize();
        
        
        try {
            String currChr = null;
            System.out.println("LCJ - byCHrom: size of infiles: " + infiles.size());
            for (Path filepath: infiles) {
                String inputFile = filepath.toString();
                int[] refChromPosList= null; int lastLinePos = 0; //lastLinePos holds the first position not already written
                boolean first = true;
                int lastStart= -1; int lastEnd= -1;
                fis = new FileInputStream(inputFile);
                System.out.println("Processing file: " + inputFile);
                scanner = new Scanner(fis);  
 
                while (scanner.hasNextLine()) {
                    String line = scanner.nextLine();
                    if ((line.isEmpty() || line.startsWith("#")) && first) continue; // skip comment lines
                    if (line.isEmpty() && !first) break;
                    String[] next = line.split("\t");
                   // if (chrName.equals(next[chrCol].toUpperCase())) {continue;} // this is header line - skip it
                    if (next[chrCol].toUpperCase().equals(chrName)) {continue;} // this is header line - skip it
                    if (first) { //initialize the reference position list                       
                        currChr = next[chrCol];
                        System.out.println("...working on chr "+currChr);
                        refChromPosList= getRefChromPosList(Integer.parseInt(next[chrCol]));
                        first= false;
                    }
 
                    //When the current line changes chromosome, finish writing end of last, reinitialize for new
                    if (!currChr.equals(next[chrCol])) {
                        // Write everything to null between the last line and the end of the chromosome
                        writeNull(lastLinePos,refChromPosList.length);
                        System.out.println(linesForChr+" lines output for chr "+currChr);
                        linesForChr= 0; lastLinePos = 0;
                        //if chromosomes not ordered 1-10, shut everything down and throw exception
                        if (Integer.parseInt(currChr)!=(Integer.parseInt(next[chrCol])-1)) {
                            Shutdown();
                            String message = "Chromosomes files must be named so they sort by chromosomes in ascending order !!\nYour chromosome " +
                               currChr + " came before chromosome " + next[chrCol];
                            throw new IllegalStateException(message);
                        }
                        currChr= next[chrCol]; lastStart= -1; lastEnd= -1;
                        refChromPosList= getRefChromPosList(Integer.parseInt(next[chrCol]));
                        System.out.println("...working on chr "+currChr);
                    }
                    int position= -1; //problems if in exponential notation (sometimes R does this) so do it as a try as double and cast if problems
                    try {
                        position= (posCol>-1)?Integer.parseInt(next[posCol]):Integer.parseInt(next[startCol]);
                        if (!oneBased()) { 
                            position++; // increment position so it matches what we have stored.
                        }
                    } catch(Exception e) {
                        System.out.println(line); 
                        position= Double.valueOf(next[posCol]).intValue();
                    }
                    int currPosOnChr=  Arrays.binarySearch(refChromPosList, position);
                    if (currPosOnChr<0) { //if position doesn't exist, throw exception
                        throw new IllegalStateException("Couldn't find pos indicated in line in reference file: "+line);
                    }
                    //Only take the first one, if position/chromosome is duplicated, as long as lines are (if range, just checks start)
                    if (lastStart==currPosOnChr) {
                        System.out.println("WARNING: Already recorded value for this position. Excluded\n"+line); continue;
                    }
                    if (range.value() && inclusive && lastEnd==currPosOnChr) { //if inclusive and lastEnd position same as start of next throw exception
                        throw new IllegalStateException("Ranges are overlapping and supposed to be inclusive!!!!!. If not inclusive, denote ranges by start, end");
                    }
                    // Write everything to null between the last line and the current line
                    // LCJ - this only writes null if the start (lastLinePos) is less than the end (currPosOnChr)
                    writeNull(lastLinePos,currPosOnChr);
                    //Write out the current line for debugging
//                    if (currPosOnChr==refChromPosList.length-1 || linesForChr>currPosOnChr || (currPosOnChr>9584780 && currPosOnChr<9584800)) {
//                        //System.out.println(currPosOnChr+"\t"+(refChromPosList.length-1)+"\t"+totalLines);
//                        System.out.println(line);
//                    }
                    if (posCol>-1) {
                        WriteValues(next); //Write out values for this line (LCJ "next" is the line we're processing)
                        lastLinePos= currPosOnChr+1; lastStart= currPosOnChr;
                        totalLines++; linesForChr++; totalNonZeroLines++;
                    } else {
                        // Looking for the value from the refChromPosList array for this chromosome.
                        // "next" holds the line we are processing.  This only works for SNP-only table.
                        // For full genome, it is a little simpler in that we go from start to end of range
                        // as it covers ALL positions.
                        int endColValue = Integer.parseInt(next[endCol]);
                        if (!oneBased()) endColValue++;
                        int end= Arrays.binarySearch(refChromPosList, endColValue);
                        if (end<0) {
                            throw new IllegalStateException("Couldn't find end pos for line: "+line);
                        }
                        end= (inclusive?end+1:end);
                        for (int i = currPosOnChr; i < end; i++) {
                            WriteValues(next); //Write out values for this line
                            totalLines++; linesForChr++; totalNonZeroLines++;
                        }
                        lastLinePos= end; //end, because even if inclusive has been transformed to exclusive
                        lastStart= currPosOnChr; //last is just for making sure that lines aren't duplicated, and checks start of range
                        lastEnd= end-1;
                    }
                } // lcj - end while loop for processing each line
                // Write everything to null between the last line and the end of the chromosome when it hits the end of the file
                writeNull(lastLinePos,refChromPosList.length);
                System.out.println(linesForChr+" lines output for chr "+currChr);
            } // lcj - end for loop processing each file
 
            Shutdown(); // close out writers
        } catch (Exception e) {
            e.printStackTrace();
            Shutdown();
            return null;
        }
        System.out.println("Process took " + (System.nanoTime() - time)/1e9 + " seconds.");
        System.out.println("LCJ - wrote the big file with totalReads: " + totalLines 
                + " totalZeroLines written: " + totalZeroLines + " totalNonZeroLines written: " + totalNonZeroLines);
        return null;
    } 
    
    private void writeNull(int start, int end) {
        try {
        for (int posIdx=start;posIdx < end; posIdx++) {
            if (realHM!=null) {for (Integer i:realWriters.keySet()) {
                realWriters.get(i).writeFloat(missToZero.value()?0:-Float.MAX_VALUE);
            }}
            if (log10HM!=null) {for (Integer i:log10Writers.keySet()) {
                    log10Writers.get(i).writeFloat(-Float.MAX_VALUE);
            }}
            if (shortHM!=null) {for (Integer i:shortWriters.keySet()) {
                shortWriters.get(i).writeShort(Short.MIN_VALUE);
            }}
            if (longHM!=null) {for (Integer i:longWriters.keySet()) {
                longWriters.get(i).writeLong(Long.MIN_VALUE);
            }}
            // lcj - you ALWAYS need the null, this worked in TE_LTRSuperFamily ... for MIchelle
            //String fam = sampLine.substring(tabPos[12] + 1, tabPos[13]) + "\n";
            //writerFam.writeChars(fam.getBytes()); 
            if (charHM!=null) {
                for (Integer i:charWriters.keySet()) {
                  charWriters.get(i).writeBytes("\n");
                }
            }
          
            if (intHM!=null) {for (Integer i:intWriters.keySet()) {
                intWriters.get(i).writeInt(Integer.MIN_VALUE);
            }}
            if (byteHM!=null) {for (Integer i:byteWriters.keySet()) {
                byteWriters.get(i).writeByte(Byte.MIN_VALUE);
            }}
            totalLines++; linesForChr++; totalZeroLines++;
        }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Problem writing nulls between "+start+" and "+end);
        }
    }
    
    private void WriteValues(String[] next) {
        try {
            if (realHM!=null) {for (Integer i:realWriters.keySet()) {
                double val= Double.MIN_VALUE; float v;
                try {
                    val = Double.parseDouble(next[i.intValue()]);
                } catch (NumberFormatException nfe) {val = Double.MIN_VALUE;}
                if (val==Double.MIN_VALUE) v= missToZero.value()?0:-Float.MAX_VALUE;
                else if (negToZero.value() && val<=0) v= 0;
                else v= (float) val;
                realWriters.get(i).writeFloat(v);
            }}
            if (log10HM!=null) {for (Integer i:log10Writers.keySet()) {
                double val= Double.MIN_VALUE;
                try {
                    val = Double.parseDouble(next[i.intValue()]);
                } catch (NumberFormatException nfe) {
                    val = Double.MIN_VALUE;
                }
                log10Writers.get(i).writeFloat(val==Double.MIN_VALUE?-Float.MAX_VALUE:(float)(-Math.log10(val)));
            }}
            if (shortHM!=null) {for (Integer i:shortWriters.keySet()) {
                short val= Short.MIN_VALUE;
                try {
                    val = Short.parseShort(next[i.intValue()]);
                } catch (NumberFormatException nfe) {
                    val = Short.MIN_VALUE;
                }
                shortWriters.get(i).writeShort(val);
            }}
            if (longHM!=null) {for (Integer i:longWriters.keySet()) {
                long val= Long.MIN_VALUE;
                try {val = Long.parseLong(next[i.intValue()]);
                } catch (NumberFormatException nfe) {val = Long.MIN_VALUE;}
                longWriters.get(i).writeLong(val);
            }}
            // lcj - always want \n !!
            if (charHM!=null) {
                for (Integer i:charWriters.keySet()) {
                  String charData = next[i.intValue()] + "\n";
                  charWriters.get(i).writeChars(charData);
                }
            }
            if (intHM!=null) {for (Integer i:intWriters.keySet()) {
                int val= Integer.MIN_VALUE;
                try {val = Integer.parseInt(next[i.intValue()]);
                } catch (NumberFormatException nfe) {val = Integer.MIN_VALUE;}
                intWriters.get(i).writeInt(val);
            }}
            if (byteHM!=null) {
                for (Integer i:byteWriters.keySet()) {
                  int val= Byte.MIN_VALUE;
                  try {
                    val = Integer.parseInt(next[i.intValue()]);                   
                  } catch (NumberFormatException nfe) {
                    val = Byte.MIN_VALUE;
                  }
                  byteWriters.get(i).writeByte((byte)val);
                }
             }
                
        } catch (Exception e) {
            for (String i:next) {System.out.println(i);}
            e.printStackTrace();
        }
    }
    
    public int[] getRefChromPosList(int chrom) {
        String refPosFile = refDir.value() + "/SNPPos_chrom" + chrom + ".txt";
        int[] refChromPosList;
        System.out.println("LCJ - GWAS ... createData - creating refChromPosList for chrom " + chrom);
        try {
            BufferedReader refPosbr = Utils.getBufferedReader(refPosFile);
            String refPosLine = null;
            ArrayList<Integer> refChromPos = new ArrayList<Integer>();
            while ((refPosLine = refPosbr.readLine()) != null) {                       
                refChromPos.add(Integer.parseInt(refPosLine));
            }  
            
            // add to primitive list - will be quicker to work
            refChromPosList = new int[refChromPos.size()];
            Iterator<Integer> iterator = refChromPos.iterator();
            for (int idx = 0; idx < refChromPosList.length; idx++)
            {
                refChromPosList[idx] = iterator.next().intValue();
            }
            
        } catch (Exception exc) {
            System.out.println("LCJ - no data for chrom " + chrom + " continuing ...");
            return null;
        }
        if (refChromPosList.length == 0) {
            System.out.println("LCJ - NO BYTES found for chrom " + chrom);
            return null;
        }
        return refChromPosList;
    }
    
    private void FindColumns(String filename) {
        // When broken by chromosome, column names shoudl be the
        // same for all files.  Send in any chromosome file and use
        // it to set the columns names
        String[] next= null;
        try { 
            FileInputStream fis = new FileInputStream(filename);
            Scanner scanner = new Scanner(fis); boolean first= true;
            while (scanner.hasNextLine()) {
                if (!first) break;
                String line = scanner.nextLine();
                if (line.isEmpty() || line.startsWith("#")) {continue;}
                next = line.split("\t");
                if (first) {
                    first= false;
                }
            }
            scanner.close();
            fis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        //get index of chromosome
        int search= -1;
        for (int col = 0; col < next.length; col++) {
            if (next[col].toLowerCase().startsWith("chr")) {search= col; break;}
        }
        if (search<0) {throw new IllegalStateException("Couldn't find chromosome in first non-empty line of input file");
        } else {
            chrCol= search;
            //chrName= next[search]; // now defaulting this to "CHR"
            System.out.println("Chr found in column "+search);
        }
        //Get position
        search= -1;
        if (range.value()) { //if range, looks for first last (inclusive) or start end (exclusive)
            for (int col = 0; col < next.length; col++) {
                if (next[col].equalsIgnoreCase("start")) {
                    search= col; break;
                }
            }
            if (search<0) {
                inclusive= true;
                for (int col = 0; col < next.length; col++) {
                    if (next[col].equalsIgnoreCase("first")) {
                        search= col; break;
                    }
                }
                if (search<0) throw new IllegalStateException("Couldn't find first or start in first line of input file to indicate range");
                else startCol= search;
            } else startCol= search;
            search= -1;
            // look for "last" or "end" for end of range
            for (int col = 0; col < next.length; col++) {
                if (next[col].equalsIgnoreCase(inclusive?"last":"end")) {
                    search= col; break;
                }
            }
            if (search<0) {
                if (search<0) throw new IllegalStateException("Couldn't find "+(inclusive?"last":"end")+" in first line of input file to indicate range");
            } else endCol= search;
            System.out.println("Using columns "+startCol+" and "+endCol+" to indicate ends of ranges, "+(inclusive?"inclusive":"exclusive"));
        } else {
            for (int col = 0; col < next.length; col++) {
                if (next[col].toLowerCase().startsWith("pos")) {
                    search= col; 
                    break;
                }
            }
            if (search<0) {System.out.println("Cannot find column for pos");
            } else {posCol= search;System.out.println("Pos found in column "+search);}
            if (chrCol<0 || posCol<0) throw new IllegalStateException("Couldn't find chromosome or position in first line of input file");
        }
        int ncols= 0;
        //Get indices of column names to keep and transform with -log10, if not null. Don't use binary search, because not sorted and won't work
        if (colsLog10.value()!=null && !colsLog10.value().equalsIgnoreCase("null") && !colsLog10.value().isEmpty()) {
            log10HM= new HashMap<>(); search= -1;
            for (String s:colsLog10.value().split(",")) {
                search= -1;
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {
                    log10HM.put(search, s);
                    log10Writers= new HashMap<>();
                    ncols++;
                    System.out.println("Found "+s+" in column "+search+" as -log10 transformed float");
                }
            }
        }
        //Get indices of column names to keep as float
        if (colsFloat.value()!=null && !colsFloat.value().equalsIgnoreCase("null") &&!colsFloat.value().isEmpty()) {
        realHM= new HashMap<>(); realWriters= new HashMap<>(); search= -1;
            for (String s:colsFloat.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {realHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as real (float)");}
            }
        }
        //Get indices of column names to keep as int
        if (colsInt.value()!=null && !colsInt.value().equalsIgnoreCase("null") && !colsInt.value().isEmpty()) {
        intHM= new HashMap<>(); intWriters= new HashMap<>(); search= -1;
            for (String s:colsInt.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {intHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as int");}
            }
        }
        //Get indices of column names to keep as short
        if (colsShort.value()!=null && !colsShort.value().equalsIgnoreCase("null") && !colsShort.value().isEmpty()) {
        shortHM= new HashMap<>(); shortWriters= new HashMap<>(); search= -1;
            for (String s:colsShort.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {shortHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as short");}
            }
        }
        //Get indices of column names to keep as long
        if (colsLong.value()!=null && !colsLong.value().equalsIgnoreCase("null") && !colsLong.value().isEmpty()) {
        longHM= new HashMap<>(); longWriters= new HashMap<>(); search= -1;
            for (String s:colsLong.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {search= col; break;}
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {longHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as long");}
            }
        }
        //Get indices of column names to keep as char
        if (colsChar.value()!=null && !colsChar.value().equalsIgnoreCase("null") && !colsChar.value().isEmpty()) {
        charHM= new HashMap<>(); charWriters= new HashMap<>(); search= -1;
            for (String s:colsChar.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {
                        search= col; break;
                    }
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {charHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as char");}
            }
        }  
        //Get indices of column names to keep as char
        if (colsByte.value()!=null && !colsByte.value().equalsIgnoreCase("null") && !colsByte.value().isEmpty()) {
        byteHM= new HashMap<>(); byteWriters= new HashMap<>(); search= -1;
            for (String s:colsByte.value().split(",")) {
                for (int col = 0; col < next.length; col++) {
                    if (s.equalsIgnoreCase(next[col])) {
                        search= col; break;
                    }
                }
                if (search<0) {System.out.println("Cannot find column "+s);
                } else {byteHM.put(search, s);ncols++;System.out.println("Found "+s+" in column "+search+" as byte");}
            }
        }  
         if (ncols<1) throw new IllegalStateException("No valid columns to read in!");   
    }
    
    private void Initialize() {
        // This method creates the .sql files for creating the table and copying into the table.
        // These values can be appended to the larger table.  Table name must be filled in by
        // user, but the skeleton of the sql command is there to feed into monetdb.
        totalLines = 0;
        totalZeroLines = 0;
        totalNonZeroLines = 0;
        linesForChr = 0;
        try {
            DataOutputStream outCopy = Utils.getDataOutputStream(outBase.value()+"_copy.sql", 1040);
            DataOutputStream outCreate = Utils.getDataOutputStream(outBase.value()+"_create.sql", 1040);
            outCopy.writeBytes("COPY binary into ... from ('./chr_int.bin','./pos_int.bin',"); boolean first= true;
            outCreate.writeBytes("CREATE TABLE ... (chr int,pos int,"); String out= new File(outBase.value()).getName(); String outDir= new File(outBase.value()).getParent();
            if (realHM!=null) {
                for (Integer i:realHM.keySet()) {
                    String outFile = out + "_"+realHM.get(i) + "_real.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'"); first= false;
                    outCreate.writeBytes("./"+out + "_"+realHM.get(i)+" real"); 
                    realWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (intHM!=null) {
                for (Integer i:intHM.keySet()) {
                    String outFile = out + "_"+ intHM.get(i) + "_int.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'"); 
                    outCreate.writeBytes("./"+out + "_"+intHM.get(i) + " int");
                    intWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (byteHM!=null) { // tinyint is a byte
                for (Integer i:byteHM.keySet()) {
                    String outFile = out + "_"+ byteHM.get(i) + "_byte.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'"); 
                    outCreate.writeBytes("./"+out + "_"+byteHM.get(i) + " tinyint"); // monetdb doesn't have "byte", tinyint instead
                    byteWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (log10HM!=null) {
                for (Integer i:log10HM.keySet()) {
                    String outFile = out + "_"+ log10HM.get(i) + "_neglog10_real.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'");
                    outCreate.writeBytes("./"+out + "_"+log10HM.get(i) + "_negLog10 real"); 
                    log10Writers.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (shortHM!=null) {
                for (Integer i:shortHM.keySet()) {
                    String outFile = out + "_"+shortHM.get(i) + "_short.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'");
                    outCreate.writeBytes("./"+out + "_"+shortHM.get(i)+" smallint"); 
                    shortWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (longHM!=null) {
                for (Integer i:longHM.keySet()) {
                    String outFile = out + "_"+longHM.get(i) + "_long.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'");
                    outCreate.writeBytes("./"+out + "_"+longHM.get(i)+" bigint"); 
                    longWriters.put(i,new LittleEndianDataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
            if (charHM!=null) {
                for (Integer i:charHM.keySet()) {
                    String outFile = out + "_"+charHM.get(i) + "_char.bin";
                    if (!first) {outCopy.writeBytes(",");outCreate.writeBytes(",");}
                    first= false;
                    outCopy.writeBytes("'"+outFile+"'");
                    outCreate.writeBytes("./"+out + "_"+charHM.get(i)+" varchar(200)"); 
                    charWriters.put(i,new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outDir+"/"+outFile))));
                }
            }
        outCopy.writeBytes(");"); outCopy.close();
        outCreate.writeBytes(");"); outCreate.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private void Shutdown() {
        try {
            scanner.close();
            fis.close();
            if (realHM!=null) {for (Integer i:realWriters.keySet()) {realWriters.get(i).close();}}
            if (intHM!=null) {for (Integer i:intWriters.keySet()) {intWriters.get(i).close();}}
            if (shortHM!=null) {for (Integer i:shortWriters.keySet()) {shortWriters.get(i).close();}}
            if (longHM!=null) {for (Integer i:longWriters.keySet()) {longWriters.get(i).close();}}
            if (charHM!=null) {for (Integer i:charWriters.keySet()) {charWriters.get(i).close();}}
            if (log10HM!=null) {for (Integer i:log10Writers.keySet()) {log10Writers.get(i).close();}}
            if (byteHM!=null) {for (Integer i:byteWriters.keySet()) {byteWriters.get(i).close();}}

        } catch (Exception e) {
            System.out.println("Problem with shutdown");
            e.printStackTrace();
        }
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
         GeneratePluginCode.generate(ColumnsToBinarySNPOnlyTablePlugin.class);
     }
     /**
      * Convenience method to run plugin with one return object.
      */
     // TODO: Replace <Type> with specific type.
//     public <Type> runPlugin(DataSet input) {
//         return (<Type>) performFunction(input).getData(0).getData();
//     }

     /**
      * File or a Directory containing Tab-delimited files with data
      * to add to the database. Files must be named chr01.txt,
      * chr02.txt etc!
      *
      * @return Input File 
      */
     public String inputFile() {
         System.out.println("LCJ - inputFile call, will return: " + inputFile.value());
         return inputFile.value();
     }

     /**
      * Set Input Directory. Directory containing Tab-delimited
      * files with data to add to the database. Files must
      * be named chr01.txt, chr02.txt etc!
      *
      * @param value Input Directory or File
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin inputFile(String value) {
         inputFile = new PluginParameter<>(inputFile, value);
         return this;
     }

     /**
      * Directory holding files that contain all the SNP physical
      * positions, separated by chr
      *
      * @return refDir
      */
     public String refDir() {
         return refDir.value();
     }

     /**
      * Set refDir. Directory holding files that contain all
      * the SNP physical positions, separated by chr
      *
      * @param value refDir
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin refDir(String value) {
         refDir = new PluginParameter<>(refDir, value);
         return this;
     }

     /**
      * Output directory and base filename to hold the binary
      * files. Will make directory if neccesary
      *
      * @return outBase
      */
     public String outBase() {
         return outBase.value();
     }

     /**
      * Set outBase. Output directory and base filename to
      * hold the binary files. Will make directory if neccesary
      *
      * @param value outBase
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin outBase(String value) {
         outBase = new PluginParameter<>(outBase, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate real
      * binaries for
      *
      * @return Columns keep as Real (Float)
      */
     public String colsFloat() {
         return colsFloat.value();
     }

     /**
      * Set Columns keep as Real (Float). Comma separated list
      * of column names to generate real binaries for
      *
      * @param value Columns keep as Real (Float)
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colsFloat(String value) {
         colsFloat = new PluginParameter<>(colsFloat, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate int
      *  binaries for
      *
      * @return Columns keep as Int
      */
     public String colsInt() {
         return colsInt.value();
     }

     /**
      * Set Columns keep as Int. Comma separated list of column
      * names to generate int  binaries for
      *
      * @param value Columns keep as Int
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colsInt(String value) {
         colsInt = new PluginParameter<>(colsInt, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate short
      * binaries for
      *
      * @return Columns keep as Short
      */
     public String colsShort() {
         return colsShort.value();
     }

     /**
      * Set Columns keep as Short. Comma separated list of
      * column names to generate short binaries for
      *
      * @param value Columns keep as Short
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colsShort(String value) {
         colsShort = new PluginParameter<>(colsShort, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate long
      * binaries for
      *
      * @return Columns keep as Long
      */
     public String colsLong() {
         return colsLong.value();
     }

     /**
      * Set Columns keep as Long. Comma separated list of column
      * names to generate long binaries for
      *
      * @param value Columns keep as Long
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colsLong(String value) {
         colsLong = new PluginParameter<>(colsLong, value);
         return this;
     }

     /**
      * Comma separated list of column names to generate byte
      * binaries for
      *
      * @return Columns keep as Byte
      */
     public String colsByte() {
         return colsByte.value();
     }

     /**
      * Set Columns keep as Byte. Comma separated list of column
      * names to generate byte binaries for
      *
      * @param value Columns keep as bhte
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colByte(String value) {
         colsByte = new PluginParameter<>(colsByte, value);
         return this;
     }
     
     /**
      * Comma separated list of column names to generate char
      * binaries for
      *
      * @return Columns keep as Char
      */
     public String colsChar() {
         return colsChar.value();
     }

     /**
      * Set Columns keep as Char. Comma separated list of column
      * names to generate char binaries for
      *
      * @param value Columns keep as Char
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colsChar(String value) {
         colsChar = new PluginParameter<>(colsChar, value);
         return this;
     }

     /**
      * Comma separated list of column names to first transform
      * using -log10 then generate binaries for
      *
      * @return Columns to Keep and transform -log10
      */
     public String colsLog10() {
         return colsLog10.value();
     }

     /**
      * Set Columns to Keep and transform -log10. Comma separated
      * list of column names to first transform using -log10
      * then generate binaries for
      *
      * @param value Columns to Keep and transform -log10
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin colsLog10(String value) {
         colsLog10 = new PluginParameter<>(colsLog10, value);
         return this;
     }

     /**
      * Columns for range data. If true, will look for 'start'
      * and 'end' (inclusive exclusive) or 'first' 'last' (inclusive
      * inclusive) instead of 'Pos' 
      *
      * @return Range information?
      */
     public Boolean range() {
         return range.value();
     }

     /**
      * Set Range information?. Columns for range data. If
      * true, will look for 'start' and 'end' (inclusive exclusive)
      * or 'first' 'last' (inclusive inclusive) instead of
      * 'Pos' 
      *
      * @param value Range information?
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin range(Boolean value) {
         range = new PluginParameter<>(range, value);
         return this;
     }

     /**
      * Will set negative float denoted column values to zero
      * (otherwise they are set to Float.MIN
      *
      * @return Negative floats to zero?
      */
     public Boolean negToZero() {
         return negToZero.value();
     }

     /**
      * Set Negative floats to zero?. Will set negative float
      * denoted column values to zero (otherwise they are set
      * to Float.MIN
      *
      * @param value Negative floats to zero?
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin negToZero(Boolean value) {
         negToZero = new PluginParameter<>(negToZero, value);
         return this;
     }

     /**
      * Will set missing float denoted column values to zero
      * (otherwise they are set to Float.MIN
      *
      * @return Missing float values to zero?
      */
     public Boolean missToZero() {
         return missToZero.value();
     }

     /**
      * Set Missing float values to zero?. Will set missing
      * float denoted column values to zero (otherwise they
      * are set to Float.MIN
      *
      * @param value Missing float values to zero?
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin missToZero(Boolean value) {
         missToZero = new PluginParameter<>(missToZero, value);
         return this;
     }
     /**
      * Will assume positions and ranges are 1-based
      * unless otherwise set.
      * 
      * @return positions are 1-based?
      */
     public Boolean oneBased() {
         return oneBased.value();
     }

     /**
      * Set positions to 1-based unless this variable
      * is "false".  
      *
      * @param value 1-based positions?
      *
      * @return this plugin
      */
     public ColumnsToBinarySNPOnlyTablePlugin oneBased(Boolean value) {
         oneBased = new PluginParameter<>(oneBased, value);
         return this;
     }
}
