/*
 *  AdjustPhasingPlugin
 * 
 *  Created on Feb 17, 2017
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
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

    private static final int NUM_VCF_HEADER_COLUMNS = 9;

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

    private PluginParameter<String> myOutputVCFFile = new PluginParameter.Builder<>("outputVCFFile", null, String.class)
            .description("Output VCF file")
            .outFile()
            .guiName("Output VCF File")
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

            Map<String, ProcessHapcut> processedHapcutFiles = new HashMap<>();
            for (Future<ProcessHapcut> future : futures) {
                ProcessHapcut processed = future.get();
                processedHapcutFiles.put(processed.taxaName(), processed);
                myLogger.info("finished processing: " + processed.filename() + " taxa: " + processed.taxaName());
            }

            try (BufferedReader reader = Utils.getBufferedReader(inputVCFFile());
                    BufferedWriter writer = Utils.getBufferedWriter(outputVCFFile())) {

                String line = reader.readLine();
                while (line != null && line.startsWith("##")) {
                    writer.write(line);
                    line = reader.readLine();
                }

                if (line == null || !line.startsWith("#CHROM")) {
                    throw new IllegalArgumentException("AdjustPhasingPlugin: processData: First line after ## lines should be header #CHROM...");
                }

                writer.write(line);
                writer.write('\n');
                String[] tokens = line.split("\t");
                Switch[] switches = new Switch[tokens.length - NUM_VCF_HEADER_COLUMNS];
                for (int i = NUM_VCF_HEADER_COLUMNS; i < tokens.length; i++) {
                    ProcessHapcut current = processedHapcutFiles.remove(tokens[i]);
                    if (current == null) {
                        switches[i - NUM_VCF_HEADER_COLUMNS] = FALSE_SWITCH;
                    } else {
                        switches[i - NUM_VCF_HEADER_COLUMNS] = current.result();
                    }
                }

                BlockingQueue<Future<ProcessLines>> queue = new LinkedBlockingQueue<>();

                Future<?> readFuture = pool.submit(new ReadLines(reader, queue, switches, pool));

                Future<?> writeFuture = pool.submit(new WriteLines(writer, queue));

                readFuture.get();

                writeFuture.get();

            } catch (Exception ex) {
                myLogger.debug(ex.getMessage(), ex);
                throw new IllegalStateException("AdjustPhasingPlugin: processData: problem converting file: " + inputVCFFile() + " to file: " + outputVCFFile() + " error: " + ex.getMessage());
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

    /**
     * Output VCF file
     *
     * @return Output VCF File
     */
    public String outputVCFFile() {
        return myOutputVCFFile.value();
    }

    /**
     * Set Output VCF File. Output VCF file
     *
     * @param value Output VCF File
     *
     * @return this plugin
     */
    public AdjustPhasingPlugin outputVCFFile(String value) {
        myOutputVCFFile = new PluginParameter<>(myOutputVCFFile, value);
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

        public String filename() {
            return myFilename;
        }

        public String taxaName() {
            return myTaxaName;
        }

        public Switch result() {
            if (myResult.cardinality() == 0) {
                return FALSE_SWITCH;
            } else {
                return new SwitchBitSet(myResult);
            }
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

    private final Switch FALSE_SWITCH = new Switch();

    private class Switch {

        public boolean alleles(int site) {
            return false;
        }
    }

    private class SwitchBitSet extends Switch {

        private final BitSet myBitSet;

        public SwitchBitSet(BitSet bitSet) {
            myBitSet = bitSet;
        }

        public boolean alleles(int site) {
            return myBitSet.fastGet(site);
        }
    }

    private class ProcessLines implements Callable<ProcessLines> {

        private final int myStartSite;
        private final String[] myLines;
        private final Switch[] mySwitches;
        private final int myNumTaxa;
        private final StringBuilder myBuilder = new StringBuilder();
        private final boolean myIsFinalProcessLines;

        public ProcessLines(int startSite, String[] lines, Switch[] switches) {
            myStartSite = startSite;
            myLines = lines;
            mySwitches = switches;
            myNumTaxa = switches.length;
            myIsFinalProcessLines = false;
        }

        public ProcessLines() {
            myStartSite = 0;
            myLines = new String[0];
            mySwitches = null;
            myNumTaxa = 0;
            myIsFinalProcessLines = true;
        }

        @Override
        public ProcessLines call() throws Exception {

            int currentSite = myStartSite;
            for (String current : myLines) {
                char[] temp = current.toCharArray();
                int numChars = temp.length;
                int copyToNthTab = NUM_VCF_HEADER_COLUMNS;
                int numTabsFound = 0;
                int startIndex = 0;
                int charIndex = 0;
                for (int t = 0; t < myNumTaxa; t++) {
                    if (mySwitches[t].alleles(currentSite)) {
                        while (true) {
                            if (temp[charIndex] == '\t') {
                                numTabsFound++;
                                if (numTabsFound == copyToNthTab) {
                                    myBuilder.append(temp, startIndex, charIndex - startIndex + 1);
                                    charIndex++;
                                    char first = temp[charIndex];
                                    charIndex++;
                                    if (temp[charIndex] != '|') {
                                        throw new IllegalStateException("AdjustPhasingPlugin: ProcessLines: | char expected.");
                                    }
                                    charIndex++;
                                    myBuilder.append(temp[charIndex]);
                                    myBuilder.append('|');
                                    myBuilder.append(first);
                                    charIndex++;
                                    while (charIndex < numChars && temp[charIndex] != '\t') {
                                        myBuilder.append(temp[charIndex]);
                                        charIndex++;
                                    }
                                    if (charIndex < numChars) {
                                        myBuilder.append('\t');
                                    } else {
                                        myBuilder.append('\n');
                                    }
                                    startIndex = charIndex + 1;
                                    break;
                                }
                            }
                            charIndex++;
                        }
                        copyToNthTab++;
                    } else {
                        copyToNthTab++;
                    }
                }

                if (startIndex < numChars) {
                    myBuilder.append(temp, startIndex, numChars - startIndex);
                    myBuilder.append('\n');
                }

                currentSite++;
            }

            return this;

        }

    }

    private static final int NUM_LINES_PER_BLOCK = 100;

    private class ReadLines implements Runnable {

        private final BufferedReader myReader;
        private final BlockingQueue<Future<ProcessLines>> myQueue;
        private final Switch[] mySwitches;
        private final ExecutorService myPool;

        public ReadLines(BufferedReader reader, BlockingQueue<Future<ProcessLines>> queue, Switch[] switches, ExecutorService pool) {
            myReader = reader;
            myQueue = queue;
            mySwitches = switches;
            myPool = pool;
        }

        @Override
        public void run() {

            try {

                List<String> temp = new ArrayList<>();
                String line = myReader.readLine();
                int count = 0;
                int startSite = 0;
                while (line != null) {
                    temp.add(line);
                    count++;
                    if (count == NUM_LINES_PER_BLOCK) {
                        String[] lines = new String[NUM_LINES_PER_BLOCK];
                        myQueue.add(myPool.submit(new ProcessLines(startSite, temp.toArray(lines), mySwitches)));
                        temp.clear();
                        count = 0;
                        startSite += NUM_LINES_PER_BLOCK;
                    }
                    line = myReader.readLine();
                }

                if (!temp.isEmpty()) {
                    String[] lines = new String[temp.size()];
                    myQueue.add(myPool.submit(new ProcessLines(startSite, temp.toArray(lines), mySwitches)));
                }

                myQueue.add(myPool.submit(new ProcessLines()));

            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: ReadLines: problem reading file: " + inputVCFFile());
            }
        }

    }

    private class WriteLines implements Runnable {

        private final BufferedWriter myWriter;
        private final BlockingQueue<Future<ProcessLines>> myQueue;

        public WriteLines(BufferedWriter writer, BlockingQueue<Future<ProcessLines>> queue) {
            myWriter = writer;
            myQueue = queue;
        }

        @Override
        public void run() {

            try {
                Future<ProcessLines> future = myQueue.take();
                ProcessLines processed = future.get();
                while (!processed.myIsFinalProcessLines) {
                    myWriter.write(processed.myBuilder.toString());
                    future = myQueue.take();
                    processed = future.get();
                }
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("AdjustPhasingPlugin: WriteLines: problem writing file: " + outputVCFFile());
            }
        }

    }

}
