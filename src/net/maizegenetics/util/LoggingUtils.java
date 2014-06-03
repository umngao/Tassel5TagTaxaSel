/*
 *  LoggingUtils
 */
package net.maizegenetics.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

/**
 *
 * @author Terry Casstevens
 */
public class LoggingUtils {

    private static final Logger myLogger = Logger.getLogger(LoggingUtils.class);
    private static final PrintStream myOriginalOutputStream = System.out;
    private static final PrintStream myOriginalErrStream = System.err;

    private static PrintStream myPrintStream;

    private LoggingUtils() {
        // Utility Class
    }

    public static void setupLogging() {
        System.setOut(myOriginalOutputStream);
        System.setErr(myOriginalErrStream);
        sendLog4jToStdout();
    }

    public static void setupLogging(PrintStream stream) {
        System.setOut(stream);
        System.setErr(stream);
        sendLog4jToStdout();
    }

    public static void setupLogfile(String logFileName) throws FileNotFoundException {
        File logFile = new File(logFileName);
        myLogger.info("Log File: " + logFile.getAbsolutePath());
        myPrintStream = new PrintStream(logFile);
        System.setOut(myPrintStream);
        System.setErr(myPrintStream);
        sendLog4jToStdout();
    }

    public static void closeLogfile() {
        if (myPrintStream != null) {
            myPrintStream.close();
        }
    }

    private static void sendLog4jToStdout() {
        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "INFO, stdout");
        props.setProperty("log4j.appender.stdout", "org.apache.log4j.ConsoleAppender");
        props.setProperty("log4j.appender.stdout.Threshold", "info");
        props.setProperty("log4j.appender.stdout.layout", "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);
    }

}
