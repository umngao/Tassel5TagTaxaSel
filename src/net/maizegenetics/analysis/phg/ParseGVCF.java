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

        // first one containing header lines (starts with ##)
        public ProcessLines(List<String> headerLines) {
            myLines = Collections.unmodifiableList(headerLines);
            myStartLine = 1;
        }

        // data lines following header lines
        public ProcessLines(int lineNum, List<String> lines) {
            myLines = Collections.unmodifiableList(lines);
            myStartLine = lineNum;
        }

        // indicates last element in queue
        public ProcessLines() {
            myLines = null;
            myStartLine = -1;
        }

        @Override
        public ProcessLines call() throws Exception {
            return null;
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
                myQueue.add(myPool.submit(new ProcessLines(lineNum, temp)));
                lineNum += temp.size();

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
