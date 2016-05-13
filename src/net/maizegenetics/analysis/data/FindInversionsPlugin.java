/*
 *  FindInversionsPlugin
 * 
 *  Created on May 9, 2016
 */
package net.maizegenetics.analysis.data;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.TreeMap;
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
import net.maizegenetics.util.TableReportBuilder;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class FindInversionsPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FindInversionsPlugin.class);

    public static enum WINDOW_UNIT {
        Sites, Positions
    };

    private PluginParameter<WINDOW_UNIT> myWindowUnit = new PluginParameter.Builder<>("windowUnit", WINDOW_UNIT.Sites, WINDOW_UNIT.class)
            .description("Window Unit")
            .range(WINDOW_UNIT.values())
            .build();

    private PluginParameter<Integer> myStepSize = new PluginParameter.Builder<>("stepSize", 100, Integer.class)
            .description("Step Size")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<Integer> myWindowSize = new PluginParameter.Builder<>("windowSize", 500, Integer.class)
            .description("Window Size")
            .range(Range.atLeast(0))
            .build();

    private PluginParameter<Integer> myNumPCAs = new PluginParameter.Builder<>("numPCAs", 1, Integer.class)
            .description("")
            .range(Range.atLeast(1))
            .guiName("Number of PCAs")
            .build();

    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<>("outputFile", null, String.class)
            .description("")
            .outFile()
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
        mds.numberOfAxes(numPCAs());

        Chromosome[] chromosomes = genotypeTable.chromosomes();

        for (Chromosome c : chromosomes) {

            FilterSiteBuilderPlugin filterChr = new FilterSiteBuilderPlugin(null, false)
                    .startChr(c)
                    .endChr(c);
            DataSet genotypeDataSet = filterChr.performFunction(input);
            GenotypeTable genotypeChr = (GenotypeTable) genotypeDataSet.getData(0).getData();
            int numSites = genotypeChr.numberOfSites();

            for (int start = 0;; start += stepSize()) {

                int startSite = -1;
                int endSite = -1;
                if (windowUnit() == WINDOW_UNIT.Sites) {
                    if (start > numSites) {
                        break;
                    }
                    startSite = start;
                    endSite = Math.min(start + windowSize() - 1, numSites - 1);
                } else if (windowUnit() == WINDOW_UNIT.Positions) {
                    startSite = genotypeChr.siteOfPhysicalPosition(start, c);
                    endSite += start + windowSize();
                    if (startSite < 0) {
                        startSite = -(startSite + 1);
                        if (startSite >= numSites) {
                            startSite = numSites - 1;
                        }
                    }
                    if (startSite > numSites) {
                        break;
                    }
                    endSite = genotypeChr.siteOfPhysicalPosition(endSite, c);
                    if (endSite < 0) {
                        endSite = -endSite - 2;
                        if (endSite >= numSites) {
                            endSite = numSites - 1;
                        }
                    }
                    if (startSite > endSite) {
                        continue;
                    }
                }

                FilterSiteBuilderPlugin filter = new FilterSiteBuilderPlugin(null, false)
                        .startSite(startSite)
                        .endSite(endSite);
                DataSet filteredGenotype = filter.performFunction(genotypeDataSet);

                DataSet distanceMatrix = distance.performFunction(filteredGenotype);

                try {
                    myLogger.info("Starting MDS...");
                    DataSet mdsResults = mds.processData(distanceMatrix);
                    myLogger.info("Finsihed MDS...");
                    for (int i = 1; i <= numPCAs(); i++) {
                        scoreSinglePCA((CorePhenotype) mdsResults.getDataOfType(CorePhenotype.class).get(0).getData(),
                                new ChrPos(c, startSite, endSite, genotypeChr.chromosomalPosition(startSite), genotypeChr.chromosomalPosition(endSite)), i);
                    }
                } catch (Exception e) {
                    myLogger.debug(e.getMessage(), e);
                    myLogger.warn("Problem calculating MDS for window chr: " + c.getName() + " start site: " + startSite + " end site: " + endSite);
                }

            }

        }

        String[] columnHeaders = new String[]{"Chromosome", "Start Position", "End Position"};
        TableReportBuilder builder = TableReportBuilder.getInstance("Candidates", columnHeaders);
        try (BufferedWriter writer = Utils.getBufferedWriter(Utils.addSuffixIfNeeded(outputFile(), ".txt"))) {
            writer.write("Chromosome\tStart Site\tEnd Site\tStart Postion\tEnd Position\tNum Gap Bins 1\tPeak Bin 1\tNum Peak Bins 1\tNum Gap Bins 2\tPeak Bin 2\tNum Peak Bins 2\tNum Gap Bins 3\tPeak Bin 3\tNum Peak Bins 3\tNum Gap Bins 4\t...\n");
            for (Map.Entry<ChrPos, List<Float>> current : myResults.entrySet()) {
                ChrPos chrPos = current.getKey();
                writer.write(chrPos.myChr.getName() + "\t" + chrPos.myStartSite + "\t" + chrPos.myEndSite + "\t" + chrPos.myStartPos + "\t" + chrPos.myEndPos);
                int numBigGaps = 0;
                int numSmallRegionPeaks = 0;
                int count = 0;
                for (Float peak : current.getValue()) {
                    int resultType = count % 3;
                    if ((resultType == 0) && (peak > 20.0f)) {
                        numBigGaps++;
                    }
                    if ((resultType == 2) && (peak < 20.0f)) {
                        numSmallRegionPeaks++;
                    }
                    count++;
                    writer.write("\t" + peak);
                }
                if ((numBigGaps >= 2) && (numSmallRegionPeaks >= 3)) {
                    Object[] row = new Object[3];
                    row[0] = chrPos.myChr;
                    row[1] = chrPos.myStartPos;
                    row[2] = chrPos.myEndPos;
                    builder.addElements(row);
                }
                writer.write("\n");
            }
        } catch (Exception e) {
            throw new IllegalStateException("Problem writing file: " + outputFile());
        }

        return new DataSet(new Datum("Candidate Inversions", builder.build(), null), this);
    }

    private void scoreSinglePCA(CorePhenotype pcaResults, ChrPos chrPos, int whichPCA) {

        int rowCount = (int) pcaResults.getRowCount();
        int numBins = Math.min(100, rowCount);
        float[] pcaValues = new float[rowCount];
        for (int i = 0; i < rowCount; i++) {
            pcaValues[i] = (float) pcaResults.getRow(i)[whichPCA];
        }

        Arrays.sort(pcaValues);
        float minValue = pcaValues[0];
        float maxValue = pcaValues[rowCount - 1];
        float range = maxValue - minValue;
        float increment = range / ((float) numBins - 1.0f);

        float step = minValue;
        int[] bins = new int[numBins];
        int count = 0;
        for (int i = 0; i < numBins - 1; i++) {
            step += increment;
            while ((count < rowCount) && (pcaValues[count] < step)) {
                bins[i]++;
                count++;
            }
        }
        bins[numBins - 1] = rowCount - count;

        int threshold = Math.round((float) rowCount * 0.005f);

        List<Float> peaksGapsWidths = new ArrayList<>();
        float currentPeak = 0.0f;
        int totalWeight = 0;
        float gap = 0.0f;
        float region = 0.0f;
        boolean inRegion = false;
        for (int i = 0; i < numBins; i++) {
            if (bins[i] <= threshold) {
                if (inRegion) {
                    currentPeak /= (float) totalWeight;
                    peaksGapsWidths.add(currentPeak);
                    currentPeak = 0.0f;
                    totalWeight = 0;
                    peaksGapsWidths.add(region);
                    region = 0.0f;
                }
                gap++;
                inRegion = false;
            } else {
                if (!inRegion) {
                    peaksGapsWidths.add(gap);
                    gap = 0.0f;
                }
                region++;
                currentPeak += (float) bins[i] * (float) i;
                totalWeight += bins[i];
                inRegion = true;
            }
        }

        if (totalWeight > 0) {
            currentPeak /= (float) totalWeight;
            peaksGapsWidths.add(currentPeak);
        }

        myResults.put(chrPos, peaksGapsWidths);

    }

    private Map<ChrPos, List<Float>> myResults = new TreeMap<>();
    private Map<ChrPos, PriorityQueue<Edge>> myResult = new TreeMap<>();

    private void scoreMDS(CorePhenotype pcaResults, ChrPos chrPos) {

        myEdges.clear();
        Cluster.clearCachedClusters();

        long rowCount = pcaResults.getRowCount();
        for (long i = 0; i < rowCount; i++) {
            Object[] current = pcaResults.getRow(i);
            Taxon taxon1 = (Taxon) current[0];
            float pca1x = (float) current[1];
            float pca1y = (float) current[2];
            Cluster first = Cluster.getInstance(taxon1, pca1x, pca1y);
            for (long j = i + 1; j < rowCount; j++) {
                Object[] next = pcaResults.getRow(j);
                Taxon taxon2 = (Taxon) next[0];
                float pca2x = (float) next[1];
                float pca2y = (float) next[2];
                Cluster second = Cluster.getInstance(taxon2, pca2x, pca2y);
                Edge edge = new Edge(first, second);
                myEdges.add(edge);
            }
        }

        reduce(3);

        for (Edge edge : myEdges) {
            System.out.println(edge.myCluster1 + "\t" + edge.myCluster2 + "\t" + edge.myDistance);
        }

        myResult.put(chrPos, new PriorityQueue<>(myEdges));

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

        int weight1 = cluster1.numTaxa();
        int weight2 = cluster2.numTaxa();
        int totalWeight = weight1 + weight2;

        Map<Cluster, Edge> edges1 = new HashMap<>(cluster1.myEdges);
        Map<Cluster, Edge> edges2 = new HashMap<>(cluster2.myEdges);

        for (Map.Entry<Cluster, Edge> entry : edges1.entrySet()) {
            Edge firstEdge = entry.getValue();
            Cluster same = entry.getKey();
            if (same == cluster2) {
                removeEdge(firstEdge);
            } else {
                Edge secondEdge = edges2.get(same);
                float distance = ((firstEdge.myDistance * (float) weight1) + (secondEdge.myDistance * (float) weight2)) / (float) totalWeight;
                Edge newEdge = new Edge(newCluster, same);
                removeEdge(firstEdge);
                removeEdge(secondEdge);
                myEdges.add(newEdge);
            }
        }

        reduceClusters(numEdges);

    }

    private void removeEdge(Edge edge) {
        edge.myCluster1.myEdges.remove(edge.myCluster2);
        edge.myCluster2.myEdges.remove(edge.myCluster1);
        myEdges.remove(edge);
    }

    PriorityQueue<Edge> myEdges = new PriorityQueue<>();

    private static class Cluster {

        private static Map<Taxon, Cluster> myInstances = new HashMap<>();

        private final List<Taxon> myList = new ArrayList<>();
        private final Map<Cluster, Edge> myEdges = new HashMap<>();
        private final float myX;
        private final float myY;

        private Cluster(Taxon taxon, float x, float y) {
            myList.add(taxon);
            myX = x;
            myY = y;
        }

        private Cluster(List<Taxon> taxa, float x, float y) {
            myList.addAll(taxa);
            myX = x;
            myY = y;
        }

        public static Cluster getInstance(Taxon taxon, float x, float y) {
            Cluster result = myInstances.get(taxon);
            if (result == null) {
                result = new Cluster(taxon, x, y);
                myInstances.put(taxon, result);
            }
            return result;
        }

        public static Cluster getInstance(Cluster cluster1, Cluster cluster2) {
            List<Taxon> result = new ArrayList<>();
            result.addAll(cluster1.myList);
            result.addAll(cluster2.myList);
            float weight1 = (float) cluster1.numTaxa();
            float weight2 = (float) cluster2.numTaxa();
            float totalWeight = weight1 + weight2;
            return new Cluster(result, (cluster1.myX * weight1 + cluster2.myX * weight2) / totalWeight, (cluster1.myY * weight1 + cluster2.myY * weight2) / totalWeight);
        }

        public int numTaxa() {
            return myList.size();
        }

        public void addEdge(Edge edge) {
            if (edge.myCluster1 != this) {
                myEdges.put(edge.myCluster1, edge);
            } else {
                myEdges.put(edge.myCluster2, edge);
            }
        }

        public static void clearCachedClusters() {
            myInstances.clear();
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

        public Edge(Cluster cluster1, Cluster cluster2) {
            myCluster1 = cluster1;
            myCluster2 = cluster2;
            myDistance = calculateDistance(cluster1.myX, cluster1.myY, cluster2.myX, cluster2.myY);
            myCluster1.addEdge(this);
            myCluster2.addEdge(this);
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
    private static float calculateDistance(float x1, float y1, float x2, float y2) {
        float xSqr = (x1 - x2) * (x1 - x2);
        float ySqr = (y1 - y2) * (y1 - y2);
        return (float) Math.sqrt(xSqr + ySqr);
    }

    private class ChrPos implements Comparable<ChrPos> {

        private final Chromosome myChr;
        private final int myStartSite;
        private final int myEndSite;
        private final int myStartPos;
        private final int myEndPos;

        public ChrPos(Chromosome chr, int startSite, int endSite, int startPos, int endPos) {
            myChr = chr;
            myStartSite = startSite;
            myEndSite = endSite;
            myStartPos = startPos;
            myEndPos = endPos;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof ChrPos)) {
                return false;
            }
            ChrPos other = (ChrPos) obj;
            if ((myChr == other.myChr) && (myStartPos == other.myStartPos)) {
                return true;
            } else {
                return false;
            }
        }

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 37 * hash + myChr.getChromosomeNumber();
            hash = 37 * hash + myStartPos;
            return hash;
        }

        @Override
        public int compareTo(ChrPos o) {
            if (myChr.getChromosomeNumber() < o.myChr.getChromosomeNumber()) {
                return -1;
            } else if (myChr.getChromosomeNumber() > o.myChr.getChromosomeNumber()) {
                return 1;
            }
            if (myStartPos < o.myStartPos) {
                return -1;
            } else if (myStartPos > o.myStartPos) {
                return 1;
            } else {
                return 0;
            }
        }

    }

    /**
     * Window Unit
     *
     * @return Window Unit
     */
    public WINDOW_UNIT windowUnit() {
        return myWindowUnit.value();
    }

    /**
     * Set Window Unit. Window Unit
     *
     * @param value Window Unit
     *
     * @return this plugin
     */
    public FindInversionsPlugin windowUnit(WINDOW_UNIT value) {
        myWindowUnit = new PluginParameter<>(myWindowUnit, value);
        return this;
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

    /**
     * Num PCAs
     *
     * @return Num PCAs
     */
    public Integer numPCAs() {
        return myNumPCAs.value();
    }

    /**
     * Set Num PCAs. Num PCAs
     *
     * @param value Num PCAs
     *
     * @return this plugin
     */
    public FindInversionsPlugin numPCAs(Integer value) {
        myNumPCAs = new PluginParameter<>(myNumPCAs, value);
        return this;
    }

    /**
     * Output File
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. Output File
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public FindInversionsPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
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
