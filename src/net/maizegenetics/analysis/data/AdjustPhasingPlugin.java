/*
 *  AdjustPhasingPlugin
 * 
 *  Created on Feb 17, 2017
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class AdjustPhasingPlugin extends AbstractPlugin {
    
    private static final Logger myLogger = Logger.getLogger(AdjustPhasingPlugin.class);
    
    private PluginParameter<String> myInputVCFFile = new PluginParameter.Builder<>("inputVCFFile", null, String.class)
            .description("Input VCF file")
            .inFile()
            .guiName("Input VCF File")
            .required(true)
            .build();
    
    private PluginParameter<String> myHapcutDir = new PluginParameter.Builder<>("hapcutDir", null, String.class)
            .description("Directory containing Hapcut output files.")
            .inDir()
            .required(true)
            .build();
    
    public AdjustPhasingPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public DataSet processData(DataSet input) {
        
        int numSites = Utils.getNumberLinesNotHashOrBlank(inputVCFFile());
        
        try {
            
            ExecutorService pool = Executors.newWorkStealingPool();
            
            List<Future<ProcessHapcut>> futures = new ArrayList<>();
            
            File temp = new File(hapcutDir());
            for (File current : temp.listFiles()) {
                String name = current.getCanonicalPath();
                if (name.endsWith("_haplotypes")) {
                    ProcessHapcut process = new ProcessHapcut(name, numSites);
                    futures.add(pool.submit(process));
                }
            }
            
            int numFutures = futures.size();
            int count = 0;
            for (Future<ProcessHapcut> future : futures) {
                ProcessHapcut processed = future.get();
                System.out.println(processed.taxaName() + " finished.");
                count++;
                progress(count * 100 / numFutures, null);
            }
            
            pool.shutdown();
            
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
        }
        
        return null;
        
    }

    /**
     * Input VCF file
     *
     * @return Input VCF File
     */
    public String inputVCFFile() {
        return myInputVCFFile.value();
    }

    /**
     * Set Input VCF File. Input file
     *
     * @param value Input VCF File
     *
     * @return this plugin
     */
    public AdjustPhasingPlugin inputVCFFile(String value) {
        myInputVCFFile = new PluginParameter<>(myInputVCFFile, value);
        return this;
    }

    /**
     * Directory containing Hapcut output files.
     *
     * @return Hapcut Dir
     */
    public String hapcutDir() {
        return myHapcutDir.value();
    }

    /**
     * Set Hapcut Dir. Directory containing Hapcut output files.
     *
     * @param value Hapcut Dir
     *
     * @return this plugin
     */
    public AdjustPhasingPlugin hapcutDir(String value) {
        myHapcutDir = new PluginParameter<>(myHapcutDir, value);
        return this;
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }
    
    @Override
    public String getButtonName() {
        return "Adjust Phasing";
    }
    
    @Override
    public String getToolTipText() {
        return "Adjust Phasing";
    }
    
    private class ProcessHapcut implements Callable<ProcessHapcut> {
        
        private final String myFilename;
        private final String myTaxaName;
        private final BitSet myResult;
        private final String[] myPositions;
        private final List<String> myChromosomes = new ArrayList<>();
        private final List<Integer> myOffsets = new ArrayList<>();
        
        public ProcessHapcut(String filename, int numSites) {
            myFilename = filename;
            myTaxaName = Utils.getFilename(filename).split("_")[0];
            myResult = new OpenBitSet(numSites);
            myPositions = new String[numSites];
        }
        
        @Override
        public ProcessHapcut call() throws Exception {
            
            try (BufferedReader reader = Utils.getBufferedReader(myFilename)) {
                
                String currentChr = null;
                int highestChrIndex = -1;
                Map<String, List<Integer>> phased = new HashMap<>();
                
                String line = reader.readLine();
                if (line == null) {
                    myLogger.warn("Hapcut file: " + myFilename + " is empty.");
                } else if (!line.startsWith("BLOCK")) {
                    throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: Expected first line to begin with BLOCK in file: " + myFilename + "\n" + line);
                }
                String blockLine = line;
                
                while (line != null) {
                    
                    if (line.startsWith("BLOCK")) {
                        
                        currentChr = null;
                        highestChrIndex = -1;
                        phased.clear();
                        blockLine = line;
                        
                    } else if (line.startsWith("********")) {
                        processBlock(currentChr, highestChrIndex, phased, blockLine);
                    } else {
                        
                        String[] tokens = line.split("\t");
                        
                        if (currentChr == null) {
                            // tokens[3] is chromosome
                            currentChr = tokens[3];
                        } else if (!currentChr.equals(tokens[3])) {
                            throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: Different chr: " + tokens[3] + " within block that started with chr: " + currentChr);
                        }

                        // tokens[0] site number plus 1
                        int tempIndex = Integer.parseInt(tokens[0]) - 1;
                        if (highestChrIndex >= tempIndex) {
                            throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: index out of order: " + tempIndex);
                        } else {
                            highestChrIndex = tempIndex;
                        }

                        // tokens[4] is physical position
                        myPositions[tempIndex] = tokens[4];

                        // tokens[1] haploid A
                        // tokens[2] haploid B
                        // tokens[7] VCF genotype field
                        String key;
                        if (tokens[7].split(":")[0].equals(tokens[1] + "|" + tokens[2])) {
                            key = "SAME";
                        } else {
                            key = "DIFFERENT";
                        }
                        List<Integer> value = phased.get(key);
                        if (value == null) {
                            List<Integer> tempList = new ArrayList<>();
                            tempList.add(tempIndex);
                            phased.put(key, tempList);
                        } else {
                            value.add(tempIndex);
                        }
                        
                    }
                    
                    line = reader.readLine();
                    
                }
                
                processBlock(currentChr, highestChrIndex, phased, blockLine);
                
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: problem reading file: " + myFilename + ": " + e.getMessage());
            }
            
            return this;
            
        }
        
        private void processBlock(String currentChr, int highestChrIndex, Map<String, List<Integer>> phased, String blockLine) {
            
            int chrIndex = myChromosomes.indexOf(currentChr);
            if (chrIndex == -1) {
                myChromosomes.add(currentChr);
                myOffsets.add(highestChrIndex);
            } else if (myOffsets.get(chrIndex) < highestChrIndex) {
                myOffsets.add(chrIndex, highestChrIndex);
            }
            
            int numStates = phased.size();
            if (numStates > 2) {
                throw new IllegalStateException("AdjustPhasingPlugin: ProcessHapcut: too many states in block: " + blockLine);
            } else if (numStates == 2) {
                List<Integer> smallest = null;
                for (List<Integer> currentList : phased.values()) {
                    if (smallest == null) {
                        smallest = currentList;
                    } else if (smallest.size() > currentList.size()) {
                        smallest = currentList;
                    }
                }
                for (Integer index : smallest) {
                    myResult.fastSet(index);
                }
            }
            
        }
        
        public String taxaName() {
            return myTaxaName;
        }
        
        public BitSet result() {
            return myResult;
        }
        
        public String[] positions() {
            return myPositions;
        }
        
        public List<String> chromosomes() {
            return myChromosomes;
        }
        
        public List<Integer> offsets() {
            return myOffsets;
        }
        
    }
    
}
