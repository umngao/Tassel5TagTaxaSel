/**
 * 
 */
package net.maizegenetics.dna.map;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;

/**
 * Utilities for reading and writing Position Lists 
 * 
 * @author lcj34
 *
 */
public class PositionListIOUtils {
    private PositionListIOUtils() {
    }
    /**
     * Returns a PositionList from a tab-delimited text SNP Conserve file.  
     * The input file has 2 tab-delimited fields indicating Chromosome Number and Position
     * A header row is the first line and looks like this:
     *   #CHROM	POS
     * The remaining rows contains integer values as below:
     *   9		18234
     * @param fileName with complete path
     * @return PositionList
     */
    public static PositionList readSNPConserveFile(String fileName) {
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            PositionListBuilder plb=new PositionListBuilder();
            //parse SNP position rows
            // First value is Chromosome number, second is position, third is Strand
            String line = fileIn.readLine(); // read/skip header
            while((line=fileIn.readLine())!=null) {
                String[] tokens=line.split("\\t");
                if (tokens.length != 2) {
                	System.err.println("Error in SNP Conserve File format:" + fileName);
                	System.err.println("Expecting tab-delimited file with 2 integer values per row");
                }
                Chromosome chrom = new Chromosome(tokens[0]);
                int pos = Integer.parseInt(tokens[1]);
                Position position = new GeneralPosition.Builder(chrom,pos).build();
                plb.add(position);
            }
            return plb.build();
        } catch(Exception e) {
            System.err.println("Error in Reading SNP Conserve File:"+fileName);
            e.printStackTrace();
        }    
        return null;
    }

}
