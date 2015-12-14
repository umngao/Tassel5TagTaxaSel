package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import org.apache.commons.math3.distribution.FDistribution;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;

public class FastMultithreadedAssociationPlugin extends AbstractPlugin {
    private GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GENOTYPE_TABLE_COMPONENT[] {
            GENOTYPE_TABLE_COMPONENT.Genotype, GENOTYPE_TABLE_COMPONENT.ReferenceProbability,
            GENOTYPE_TABLE_COMPONENT.AlleleProbability };


    //plugin parameter definitions
    private PluginParameter<Double> maxp =
            new PluginParameter.Builder<>("MaxPValue", .001, Double.class)
                    .guiName("MaxPValue")
                    .description("The maximum p-value that will be output by the analysis.")
                    .build();
    private PluginParameter<GENOTYPE_TABLE_COMPONENT> myGenotypeTable =
            new PluginParameter.Builder<>("genotypeComponent", GENOTYPE_TABLE_COMPONENT.Genotype, GENOTYPE_TABLE_COMPONENT.class)
                    .genotypeTable()
                    .range(GENOTYPE_COMP)
                    .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
                    .build();
    private PluginParameter<Boolean> saveAsFile =
            new PluginParameter.Builder<>("writeToFile", false, Boolean.class)
                    .description("Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirements"
                            + "and the results will NOT be saved to the data tree. Default = false.")
                    .guiName("Write to file")
                    .build();
    private PluginParameter<String> reportFilename =
            new PluginParameter.Builder<>("outputFile", null, String.class)
                    .outFile()
                    .dependentOnParameter(saveAsFile)
                    .description("The name of the file to which these results will be saved.")
                    .guiName("Output File")
                    .build();
    private PluginParameter<Integer> maxThreads = new PluginParameter.Builder<>("maxThreads", TasselPrefs.getMaxThreads(), Integer.class)
    		.description("the maximum number of threads to be used by this plugin.")
    		.guiName("Max Threads")
    		.build();
    
    public FastMultithreadedAssociationPlugin() {
        this(null, false);
    }

    public FastMultithreadedAssociationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        // TODO finish implementing
    	
        return null;
    }

    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        return "Fast-MT Association";
    }

    @Override
    public String getToolTipText() {
        return "Multi-threaded version of Fast Association";
    }

    class SiteTester extends Thread {
        final List<double[]> orthogonalPhenotypes;
        final List<String> phenotypeNames;
        final BlockingQueue<Marker> siteQueue;
        final BlockingQueue<Object[]> outQueue;
        final double minR2;
        final int nphenotypes;
        final int numberOfObservations;
        final double errdf;
        final private FDistribution Fdist;
        
        SiteTester(List<double[]> orthogonalPhenotypes, List<String> phenotypeNames, BlockingQueue<Marker> siteQueue, BlockingQueue<Object[]> outQueue, double minRSquare, double errdf, int ntaxa) {
            this.orthogonalPhenotypes = orthogonalPhenotypes;
            this.phenotypeNames = phenotypeNames;
            this.siteQueue = siteQueue;
            this.outQueue = outQueue;
            minR2 = minRSquare;
            numberOfObservations = ntaxa;
            this.errdf = errdf;
            nphenotypes = orthogonalPhenotypes.size();
            Fdist = new FDistribution(1, errdf);
        }
        
        public void run() {
            try {
            	Marker thisMarker = siteQueue.poll(1, TimeUnit.SECONDS);
                double[] siteValues = thisMarker.values;
                while (siteValues.length > 0) {
                    double[] r2values = new double[nphenotypes];
                    
                    for (int p = 0; p < nphenotypes; p++) {
                        int sumprod = 0;
                        double[] pheno = orthogonalPhenotypes.get(p);
                        for (int t = 0; t < numberOfObservations; t++) sumprod += siteValues[t] * pheno[t];
                        r2values[p] = sumprod;
                    }
                    
                    outputResult(r2values, thisMarker.myPosition);
                    thisMarker = siteQueue.poll(1, TimeUnit.SECONDS);
                }
            } catch (InterruptedException e) {
                throw new RuntimeException("InterruptedException occurred in SiteTester thread", e);
            }
        }
        
        private void outputResult(double[] rvalues, Position pos) throws InterruptedException {
            //for r2 values >= minR2, create an output record and add it to the output queue
        	for (int p = 0; p < nphenotypes; p++) {
        		if (rvalues[p] >= minR2) {
        			Object[] result = new Object[]{ phenotypeNames.get(p), pos.getSNPID(),
                            pos.getChromosome().getName(), pos.getPosition(),
                            1, rvalues[p],
                            pvalue(rvalues[p])};
        			outQueue.put(result);
        		}
        	}
        }
        
        private double pvalue(double rvalue) {
            double F = rvalue / (1 - rvalue) * errdf;
            double p;
            try {
                p = 1 - Fdist.cumulativeProbability(F);
            } catch (Exception e) {
                p = Double.NaN;
            }
            return p;
        }

    }
    
    class Marker {
    	double[] values;
    	Position myPosition;
    	
    	Marker(double[] values, Position pos) {
    		this.values = values;
    		myPosition = pos;
    	}
    	
    }
}
