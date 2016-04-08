/*
 *  MaskGenotypePlugin
 * 
 *  Created on May 8, 2015
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.net.URL;
import javax.swing.*;
import java.util.List;
import java.util.Random;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.MaskGenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import org.apache.log4j.Logger;

/**
 * The purpose of this class is to mask (make UNKNOWN) a portion of the
 * genotypes of the given GenotypeTable. The resulting GenotypeTable will be the
 * same except some genotypes will have been changed to Unknown.
 *
 * @author Terry Casstevens
 */
public class MaskGenotypePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MaskGenotypePlugin.class);

    private PluginParameter<Double> myPercentageMasked
            = new PluginParameter.Builder<>("percentageMasked", 0.01, Double.class)
            .range(Range.closed(0.0, 1.0))
            .description("Percentage of genotypes (not already unknown) to mask.")
            .build();

    public MaskGenotypePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("MaskGenotypePlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum inputDatum = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable original = (GenotypeTable) inputDatum.getData();
        GenotypeCallTable origCalls = original.genotypeMatrix();

        MaskGenotypeTableBuilder builder = new MaskGenotypeTableBuilder(original);

        percentageMask(builder, origCalls);

        GenotypeTable result = builder.build();

        Datum genotype = new Datum(inputDatum.getName() + "_Masked", result, null);

        return new DataSet(genotype, this);

    }

    private void percentageMask(MaskGenotypeTableBuilder builder, GenotypeCallTable origCalls) {
        
        if (percentageMasked() == 0.0) {
            return;
        }

        Random random = new Random();
        int numTaxa = origCalls.numberOfTaxa();
        int numSites = origCalls.numberOfSites();
        long totalGenotypes = numTaxa * numSites;
        long numKnown = 0;
        for (int s = 0; s < numSites; s++) {
            byte[] genotypes = origCalls.genotypeForAllTaxa(s);
            for (byte current : genotypes) {
                if (current != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    numKnown++;
                }
            }
        }

        long numberMasked = 0;
        long numberToMask = (long) Math.floor((double) numKnown * percentageMasked());
        myLogger.info("Number of Genotypes Masked: " + numberToMask);

        double step = (double) totalGenotypes / (double) numberToMask * 1.8;
        long index = 0;
        while (true) {
            index += (long) (step * random.nextDouble());
            if (index >= totalGenotypes) {
                index %= totalGenotypes;
            }
            int t = (int) (index % numTaxa);
            int site = (int) (index / numTaxa);
            if ((!builder.get(t, site))
                    && (origCalls.genotype(t, site) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)) {
                builder.set(t, site);
                numberMasked++;
                if (numberMasked >= numberToMask) {
                    return;
                }
            }
        }

    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Percentage of genotypes to mask.
     *
     * @return Percentage Masked
     */
    public Double percentageMasked() {
        return myPercentageMasked.value();
    }

    /**
     * Set Percentage Masked. Percentage of genotypes to mask.
     *
     * @param value Percentage Masked
     *
     * @return this plugin
     */
    public MaskGenotypePlugin percentageMasked(Double value) {
        myPercentageMasked = new PluginParameter<>(myPercentageMasked, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Mask Genotype";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = MaskGenotypePlugin.class.getResource("/net/maizegenetics/analysis/images/mask.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Mask Genotype";
    }

}
