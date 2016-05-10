/*
 *  FindInversionsPlugin
 * 
 *  Created on May 9, 2016
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.CorePhenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

/**
 *
 * @author Terry Casstevens
 */
public class FindInversionsPlugin extends AbstractPlugin {

    private PluginParameter<Integer> myStepSize = new PluginParameter.Builder<>("stepSize", 100, Integer.class)
            .description("Step Size")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<Integer> myWindowSize = new PluginParameter.Builder<>("windowSize", 500, Integer.class)
            .description("Window Size")
            .range(Range.atLeast(0))
            .build();

    public FindInversionsPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> data = input.getDataOfType(GenotypeTable.class);
        if (data.size() != 1) {
            throw new IllegalArgumentException("FindInversions: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum data = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable) data.getData();

        PrincipalComponentsPlugin pca = new PrincipalComponentsPlugin(getParentFrame(), false);
        pca.numberOfComponents(2);

        Chromosome[] chromosomes = genotypeTable.chromosomes();

        for (Chromosome c : chromosomes) {

            FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                    .startChr(c)
                    .endChr(c);
            DataSet genotypeDataSet = filterChr.performFunction(input);
            GenotypeTable genotypeChr = (GenotypeTable) genotypeDataSet.getData(0).getData();
            int numSites = genotypeChr.numberOfSites();

            for (int startSite = 0; startSite < numSites; startSite += stepSize()) {

                int endSite = Math.min(startSite + windowSize() - 1, numSites - 1);

                FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false)
                        .startSite(startSite)
                        .endSite(endSite);
                DataSet filteredGenotype = filter.performFunction(genotypeDataSet);

                DataSet pcaResults = pca.performFunction(filteredGenotype);
                scorePCA((CorePhenotype) pcaResults.getDataOfType(CorePhenotype.class).get(0).getData());

            }

        }

        return null;
    }

    private void scorePCA(CorePhenotype pcaResults) {
    /**
     * Calculates distance between two points.
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    private float calculateDistance(float x1, float y1, float x2, float y2) {
        return 0.0f;
    }


    }

    /**
     * Step Size
     *
     * @return Step Size
     */
    public Integer stepSize() {
        return myStepSize.value();
    }

    /**
     * Set Step Size. Step Size
     *
     * @param value Step Size
     *
     * @return this plugin
     */
    public FindInversionsPlugin stepSize(Integer value) {
        myStepSize = new PluginParameter<>(myStepSize, value);
        return this;
    }

    /**
     * Window Size
     *
     * @return Window Size
     */
    public Integer windowSize() {
        return myWindowSize.value();
    }

    /**
     * Set Window Size. Window Size
     *
     * @param value Window Size
     *
     * @return this plugin
     */
    public FindInversionsPlugin windowSize(Integer value) {
        myWindowSize = new PluginParameter<>(myWindowSize, value);
        return this;
    }

    @Override
    public String getToolTipText() {
        return "Find Inversions";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Find Inversions";
    }

    @Override
    public String getCitation() {
        return "Mei W, Casstevens T. (May 2016) Third Tassel Hackathon.";
    }

}
