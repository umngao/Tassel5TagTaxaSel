package net.maizegenetics.analysis.phg;

/**
 * @author Terry Casstevens Created June 28, 2017
 */

import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;


public class ParseGVCF {

    private static final Logger myLogger = Logger.getLogger(ParseGVCF.class);

    private static final int NUM_LINES_PER_BLOCK = 10;

    private ParseGVCF() {
        // utility
    }

    public static BlockingQueue<Future<ProcessLines>> parse(String filename) {
        ExecutorService pool = Executors.newWorkStealingPool();
        BlockingQueue<Future<ProcessLines>> queue = new LinkedBlockingQueue<>();
        pool.submit(new ReadLines(filename, queue, pool));
        return queue;
    }

    public static class ProcessLines implements Callable<ProcessLines> {

        private final List<String> myLines;
        private final int myStartLine;
        private final int myNumLines;
        private final boolean myIsHeader;
        private List<GVCFLine> myProcessedLines = null;

        // first one containing header lines (starts with #)
        public ProcessLines(List<String> headerLines) {
            myLines = Collections.unmodifiableList(headerLines);
            myStartLine = 1;
            myNumLines = headerLines.size();
            myIsHeader = true;
        }

        // data lines following header lines
        public ProcessLines(int lineNum, List<String> lines) {
            myLines = Collections.unmodifiableList(lines);
            myStartLine = lineNum;
            myNumLines = lines.size();
            myIsHeader = false;
        }

        // indicates last element in queue
        public ProcessLines() {
            myLines = null;
            myStartLine = -1;
            myNumLines = 0;
            myIsHeader = false;
        }

        @Override
        public ProcessLines call() throws Exception {

            if (isFinal() || isHeader()) {
                return this;
            }

            List<GVCFLine> temp = new ArrayList<>();
            int currentLineNum = myStartLine;
            for (String current : myLines) {
                temp.add(new GVCFLine(currentLineNum, current));
                currentLineNum++;
            }
            myProcessedLines = temp;

            return this;

        }

        public List<String> lines() {
            return myLines;
        }

        public int startLine() {
            return myStartLine;
        }

        public int numLines() {
            return myNumLines;
        }

        public boolean isHeader() {
            return myIsHeader;
        }

        public boolean isFinal() {
            return myStartLine == -1;
        }

        public List<GVCFLine> processedLines() {
            return myProcessedLines;
        }

    }

    public static class GVCFLine {

        private final String myLine;
        private final int myLineNum;
        // DP value
        private int myDepth = 0;

        public GVCFLine(int lineNum, String line) {
            myLineNum = lineNum;
            myLine = line;
            parseLine();
        }

        public int lineNum() {
            return myLineNum;
        }

        @Override
        public String toString() {
            return myLine;
        }

        public int depth() {
            return myDepth;
        }

        private void parseLine() {

            if (myLine == null || myLine.startsWith("#")) {
                myLogger.error("line: " + myLine);
                throw new IllegalStateException("ParseGVCF: line shouldn't be null or start with #");
            }

            // 0        1       2       3       4       5       6       7       8       9
            // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  W22
            // 8       17639   .       C       <NON_REF>       .       .       END=17789       GT:DP:GQ:MIN_DP:PL      0:12:99:8:0,140
            // 8       17790   .       T       A,<NON_REF>     175.00  .       DP=5;MLEAC=1,0;MLEAF=1.000,0.000;RAW_MQ=18000.00        GT:AD:DP:GQ:PL:SB       1:0,5,0:5:99:205,0,205:0,0,2,3
            String[] tokens = myLine.split("\t");

            String chr = tokens[0];
            int start = Integer.parseInt(tokens[1]);
            int end = start;

            String id = tokens[2];
            String ref = tokens[3];
            String alt = tokens[4];
            String qual = tokens[5];
            String filter = tokens[6];

            String[] info = tokens[7].split(";");
            for (String current : info) {
                String[] keyValue = current.split("=");
                if (keyValue[0].equals("END")) {
                    end = Integer.parseInt(keyValue[1]);
                }
            }

            // TODO: making end exclusive
            end++;

            int length = end - start;

            if (!alt.equals("<NON_REF>")) {
                // TODO
            }

            String[] formats = tokens[8].split(":");
            String[] values = tokens[9].split(":");

            int numFormats = formats.length;
            if (numFormats != values.length) {
                throw new IllegalArgumentException("unequal formats");
            }

            boolean isHomo = true;
            boolean isHomoBasedOnGenotype = false;
            boolean isAlt = false;
            boolean isDiploid = false;
            for (int i = 0; i < numFormats; i++) {
                if (formats[i].equals("DP")) {
                    myDepth = Integer.parseInt(values[i]);
                } else if (formats[i].equals("AD")) {
                    String[] alleleDepths = values[i].split(",");
                    // TODO - compare sum AD to DP
                    if (alleleDepths.length > 1 && !alleleDepths[0].equals("0") && !alleleDepths[1].equals("0")) {
                        isHomo = false;
                    }
                } else if (formats[i].equals("GT")) {
                    String[] alleles = values[i].split("/");
                    if (alleles.length == 1) {
                        if (!values[i].equals("0")) {
                            isAlt = true;
                        }
                    } else {
                        isDiploid = true;
                        if (alleles[0].equals(alleles[1])) {
                            isHomoBasedOnGenotype = true;
                        }
                        if (!alleles[0].equals("0") || !alleles[1].equals("0")) {
                            isAlt = true;
                        }
                    }
                }
            }

            if (isDiploid && (isHomo != isHomoBasedOnGenotype)) {
                myLogger.error("Line: " + myLine);
                throw new IllegalStateException("ParseGVCF: Homozygous doesn't match based on DP and GT");
            }

            // TODO - Unnecessary
            if (isDiploid) {
                isHomo = isHomoBasedOnGenotype;
            }

            if (myDepth > 0) {
                if (isHomo) {
                    if (isAlt) {
                        // numHomozygousAlt += length;
                        if (myDepth == 1) {
                            // numHomozygousAlt1 += length;
                        } else if (myDepth == 2) {
                            // numHomozygousAlt2 += length;
                        } else {
                            // numHomozygousAlt3 += length;
                        }
                    } else {
                        // numHomozygous += length;
                    }
                } else {
                    // numHeterozygous += length;
                }
            }

        }

    }

    private static class ReadLines implements Runnable {

        private final String myFilename;
        private final BlockingQueue<Future<ProcessLines>> myQueue;
        private final ExecutorService myPool;

        public ReadLines(String filename, BlockingQueue<Future<ProcessLines>> queue, ExecutorService pool) {
            myFilename = filename;
            myQueue = queue;
            myPool = pool;
        }

        @Override
        public void run() {

            try (BufferedReader reader = Utils.getBufferedReader(myFilename)) {

                int lineNum = 1;
                List<String> temp = new ArrayList<>();
                String line = reader.readLine();

                while (line != null && line.startsWith("#")) {
                    temp.add(line);
                    line = reader.readLine();
                }
                if (temp.size() != 0) {
                    myQueue.add(myPool.submit(new ProcessLines(temp)));
                    lineNum += temp.size();
                }

                int count = 0;
                temp = new ArrayList<>();
                while (line != null) {
                    temp.add(line);
                    count++;
                    if (count == NUM_LINES_PER_BLOCK) {
                        myQueue.add(myPool.submit(new ProcessLines(lineNum, temp)));
                        temp = new ArrayList<>();
                        count = 0;
                        lineNum += NUM_LINES_PER_BLOCK;
                    }
                    line = reader.readLine();
                }

                if (!temp.isEmpty()) {
                    String[] lines = new String[temp.size()];
                    myQueue.add(myPool.submit(new ProcessLines(lineNum, temp)));
                }

                myQueue.add(myPool.submit(new ProcessLines()));

            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("ParseGVCF: ReadLines: problem reading file: " + myFilename);
            }

        }

    }

}
