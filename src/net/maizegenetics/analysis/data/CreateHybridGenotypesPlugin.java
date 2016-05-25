/*
 *  CreateHybridGenotypesPlugin
 * 
 *  Created on May 17, 2016
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.io.BufferedReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.HybridGenotypeCallTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class CreateHybridGenotypesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(CreateHybridGenotypesPlugin.class);

    private PluginParameter<String> myHybridFile = new PluginParameter.Builder<>("hybridFile", null, String.class)
            .description("Two column tab-delimited file defining parent crosses.")
            .required(true)
            .inFile()
            .build();

    private PluginParameter<String> myHybridChar = new PluginParameter.Builder<>("hybridChar", "/", String.class)
            .description("String used to combine taxa names to create hybrid name.")
            .build();

    public CreateHybridGenotypesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("CreateHybridGenotypesPlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum data = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable) data.getData();
        TaxaList origTaxa = genotypeTable.taxa();

        TaxaListBuilder taxa = new TaxaListBuilder();
        try (BufferedReader reader = Utils.getBufferedReader(hybridFile())) {
            List<Integer> firstParents = new ArrayList<>();
            List<Integer> secondParents = new ArrayList<>();
            int lineNum = 1;
            String line = reader.readLine();
            while (line != null) {
                String[] tokens = line.trim().split("\t");
                if (tokens.length != 2) {
                    throw new IllegalArgumentException("CreateHybridGenotypePlugin: processData: Must have two parents per line in file: "
                            + hybridFile() + ". Problem on line: " + lineNum);
                }
                int firstIndex = origTaxa.indexOf(tokens[0]);
                int secondIndex = origTaxa.indexOf(tokens[1]);
                if (firstIndex == -1) {
                    myLogger.warn("processData: line: " + lineNum + " taxon name: " + tokens[0] + " not in the genotype.");
                } else if (secondIndex == -1) {
                    myLogger.warn("processData: line: " + lineNum + " taxon name: " + tokens[1] + " not in the genotype.");
                } else {
                    firstParents.add(firstIndex);
                    secondParents.add(secondIndex);
                    taxa.add(new Taxon(tokens[0] + hybridChar() + tokens[1]));
                }
                line = reader.readLine();
                lineNum++;
            }

            GenotypeTable result = GenotypeTableBuilder.getInstance(new HybridGenotypeCallTable(genotypeTable.genotypeMatrix(), firstParents, secondParents), genotypeTable.positions(), taxa.build());
            return new DataSet(new Datum(data.getName() + "_Hybrids", result, null), this);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("CreateHybridGenotypePlugin: processData: problem reading hybrid file: " + hybridFile());
        }

    }

    /**
     * Hybrid File
     *
     * @return Hybrid File
     */
    public String hybridFile() {
        return myHybridFile.value();
    }

    /**
     * Set Hybrid File. Hybrid File
     *
     * @param value Hybrid File
     *
     * @return this plugin
     */
    public CreateHybridGenotypesPlugin hybridFile(String value) {
        myHybridFile = new PluginParameter<>(myHybridFile, value);
        return this;
    }

    /**
     * String used to combine taxa names to create hybrid name.
     *
     * @return Cross Char
     */
    public String hybridChar() {
        return myHybridChar.value();
    }

    /**
     * Set Cross Char. String used to combine taxa names to create hybrid name.
     *
     * @param value Cross Char
     *
     * @return this plugin
     */
    public CreateHybridGenotypesPlugin hybridChar(String value) {
        myHybridChar = new PluginParameter<>(myHybridChar, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Create Hybrid Genotypes";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = MaskGenotypePlugin.class.getResource("/net/maizegenetics/analysis/images/hybrid.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Create Hybrid Genotypes";
    }
}
