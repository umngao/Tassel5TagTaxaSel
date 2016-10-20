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

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {

            writer.write("<Phenotype>\n");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if (i != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(phenotype.attributeType(i).name());
            }
            writer.write("\n");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if (i != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(phenotype.attributeName(i));
            }
            writer.write("\n");

            TableReportUtils.saveDelimitedTableReport(phenotype, DELIMITER, writer, false);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("PhenotypeUtils: write: problem saving file: " + filename);
        }

    }

    public static void writePlink(Phenotype phenotype, String filename) {

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {

            writer.write("FID");
            writer.write(DELIMITER);
            writer.write("IID");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if ((phenotype.attributeType(i) == Phenotype.ATTRIBUTE_TYPE.data)
                        || (phenotype.attributeType(i) == Phenotype.ATTRIBUTE_TYPE.covariate)) {
                    writer.write(DELIMITER);
                    writer.write(phenotype.attributeName(i));
                }
            }
            writer.write("\n");

            int numObservations = phenotype.numberOfObservations();
            for (int i = 0; i < numObservations; i++) {
                String taxonName = phenotype.value(i, 0).toString();
                writer.write(taxonName);
                writer.write(DELIMITER);
                writer.write(taxonName);
                writer.write(DELIMITER);
                for (int j = 1; j < phenotype.numberOfAttributes(); j++) {
                    if ((phenotype.attributeType(j) == Phenotype.ATTRIBUTE_TYPE.data)
                            || (phenotype.attributeType(j) == Phenotype.ATTRIBUTE_TYPE.covariate)) {
                        writer.write(DELIMITER);
                        String value = phenotype.value(i, j).toString();
                        if (value.equalsIgnoreCase("NaN")) {
                            writer.write("NA");
                        } else {
                            writer.write(value);
                        }
                    }
                }
                writer.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("PhenotypeUtils: writePlink: problem saving file: " + filename);
        }

    }

}
