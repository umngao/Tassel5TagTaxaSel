/*
 *  PhenotypeUtils
 * 
 *  Created on Oct 27, 2014
 */
package net.maizegenetics.phenotype;

import java.io.BufferedWriter;

import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class PhenotypeUtils {

    private static final Logger myLogger = Logger.getLogger(PhenotypeUtils.class);

    private static final String DELIMITER = "\t";

    private PhenotypeUtils() {
        // utility
    }

    public static void write(Phenotype phenotype, String filename) {

        BufferedWriter writer = null;
        try {
            writer = Utils.getBufferedWriter(filename);

            Object[] colNames = phenotype.getTableColumnNames();
            colNames[0] = "<Trait>";
            for (int j = 0; j < colNames.length; j++) {
                if (j != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(colNames[j].toString());
            }
            writer.write("\n");

            TableReportUtils.saveDelimitedTableReport(phenotype, DELIMITER, writer, false);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("PhenotypeUtils: write: problem saving file: " + filename);
        } finally {
            if (writer != null) {
                try {
                    writer.close();
                } catch (Exception ex) {
                    // do nothing
                }
            }
        }
    }

}
