/*
 *  FindInversionsPlugin
 * 
 *  Created on May 9, 2016
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.distance.DistanceMatrixPlugin;
import net.maizegenetics.analysis.distance.MultiDimensionalScalingPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.CorePhenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;

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
            throw new IllegalArgumentException("FindInversionsPlugin: preProcessParameters: must input 1 GenotypeTable.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        Datum data = input.getDataOfType(GenotypeTable.class).get(0);
        GenotypeTable genotypeTable = (GenotypeTable) data.getData();

        DistanceMatrixPlugin distance = new DistanceMatrixPlugin(getParentFrame(), false);

        MultiDimensionalScalingPlugin mds = new MultiDimensionalScalingPlugin(getParentFrame(), false);
        mds.numberOfAxes(2);

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

                DataSet distanceMatrix = distance.performFunction(filteredGenotype);

                DataSet mdsResults = mds.performFunction(distanceMatrix);
                scoreMDS((CorePhenotype) mdsResults.getDataOfType(CorePhenotype.class).get(0).getData());

            }

        }

        return null;
    }

    private void scoreMDS(CorePhenotype pcaResults) {

        long rowCount = pcaResults.getRowCount();
        for (long i = 0; i < rowCount; i++) {
            Object[] current = pcaResults.getRow(i);
            Taxon taxon1 = (Taxon) current[0];
            float pca1x = (float) current[1];
            float pca1y = (float) current[2];
            Cluster first = Cluster.getInstance(taxon1);
            for (long j = i + 1; j < rowCount; j++) {
                Object[] next = pcaResults.getRow(j);
                Taxon taxon2 = (Taxon) next[0];
                float pca2x = (float) next[1];
                float pca2y = (float) next[2];
                Cluster second = Cluster.getInstance(taxon2);
                myEdges.add(new Edge(first, second, calculateDistance(pca1x, pca1y, pca2x, pca2y)));
            }
        }

        reduce(3);

    }

    private void reduce(int numClusters) {
        int numEdges = numClusters * (numClusters + 1) / 2 - numClusters;
        reduceClusters(numEdges);
    }

    private void reduceClusters(int numEdges) {

        if (myEdges.size() <= numEdges) {
            return;
        }

        Edge current = myEdges.remove();
        Cluster cluster1 = current.myCluster1;
        Cluster cluster2 = current.myCluster2;
        Cluster newCluster = Cluster.getInstance(cluster1, cluster2);

        List<Edge> edges1 = cluster1.myEdges;
        List<Edge> edges2 = cluster2.myEdges;

    }

    PriorityQueue<Edge> myEdges = new PriorityQueue<>();

    private static class Cluster {

        private static Map<Taxon, Cluster> myInstances = new HashMap<>();

        private final List<Taxon> myList = new ArrayList<>();
        private final List<Edge> myEdges = new ArrayList<>();

        private Cluster(Taxon taxon) {
            myList.add(taxon);
        }

        private Cluster(List<Taxon> taxa) {
            myList.addAll(taxa);
        }

        public static Cluster getInstance(Taxon taxon) {
            Cluster result = myInstances.get(taxon);
            if (result == null) {
                result = new Cluster(taxon);
                myInstances.put(taxon, result);
            }
            return result;
        }

        public static Cluster getInstance(Cluster cluster1, Cluster cluster2) {
            List<Taxon> result = new ArrayList<>();
            result.addAll(cluster1.myList);
            result.addAll(cluster2.myList);
            return new Cluster(result);
        }

        public int numTaxa() {
            return myList.size();
        }

        public void addEdge(Edge edge) {
            myEdges.add(edge);
        }

        @Override
        public String toString() {
            return "Cluster{" + "myList=" + myList + '}';
        }

    }

    private static class Edge implements Comparable<Edge> {

        private final Cluster myCluster1;
        private final Cluster myCluster2;
        private final float myDistance;

        public Edge(Cluster cluster1, Cluster cluster2, float distance) {
            myCluster1 = cluster1;
            myCluster2 = cluster2;
            myDistance = distance;
        }

        @Override
        public int compareTo(Edge o) {
            if (myDistance < o.myDistance) {
                return -1;
            } else if (myDistance > o.myDistance) {
                return 1;
            } else {
                return 0;
            }
        }

    }

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
        float xSqr = (float) Math.pow(x1 - x2, 2);
        float ySqr = (float) Math.pow(y1 - y2, 2);
        return (float) Math.sqrt(xSqr + ySqr);
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
