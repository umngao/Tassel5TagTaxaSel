package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

import javax.swing.ImageIcon;

import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

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
    private PluginParameter<Boolean> addOnly =
            new PluginParameter.Builder<>("addOnly", false, Boolean.class)
                    .description("Should an additive only model be fit? If true, an additive model will be fit. If false, an additive + dominance model will be fit. Default = false.")
                    .guiName("Additive Only Model")
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
    
    public FastMultithreadedAssociationPlugin() {
        // TODO Auto-generated constructor stub
    }

    public FastMultithreadedAssociationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        // TODO Auto-generated constructor stub
    }

    @Override
    public DataSet processData(DataSet input) {
        // TODO Auto-generated method stub
        return super.processData(input);
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
        final BlockingQueue<double[]> siteQueue;
        final double minR2;
        final int nphenotypes;
        
        SiteTester(List<double[]> orthogonalPhenotypes, BlockingQueue<double[]> siteQueue, BlockingQueue<Object[]> outQueue, double minRSquare) {
            this.orthogonalPhenotypes = orthogonalPhenotypes;
            this.siteQueue = siteQueue;
            minR2 = minRSquare;
            nphenotypes = orthogonalPhenotypes.size();
        }
        
        public void run() {
            try {
                double[] siteValues = siteQueue.poll(1, TimeUnit.SECONDS);
                int ntaxa = siteValues.length;
                while (siteValues.length > 0) {
                    double[] r2values = new double[nphenotypes];
                    
                    for (int p = 0; p < nphenotypes; p++) {
                        int sumprod = 0;
                        double[] pheno = orthogonalPhenotypes.get(p);
                        for (int t = 0; t < ntaxa; t++) sumprod += siteValues[t] * pheno[t];
                        r2values[p] = sumprod;
                    }
                    
                    outputResult(r2values);
                    siteValues = siteQueue.poll(1, TimeUnit.SECONDS);
                }
            } catch (InterruptedException e) {
                throw new RuntimeException("Error polling the site queue", e);
            }
        }
        
        private void outputResult(double[] r2Values) {
            //for r2 values >= minR2, create an output record and add it to the output queue
        }
    }
}
