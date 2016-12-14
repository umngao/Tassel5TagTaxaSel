/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.ImageIcon;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.log4j.Logger;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

/**
 * @author lcj34
 *
 */
public class RepGenLDAnalysisPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(RepGenLDAnalysisPlugin.class);
    
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    
    public RepGenLDAnalysisPlugin() {
        super(null, false);
    }

    public RepGenLDAnalysisPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RepGenLDAnalysisPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenLDAnalysisPlugin: postProcessParameters: Input DB not set or found");
        }
    }
    
    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
 
        try {           
            System.out.println("RepGenLDAnalysis:processData begin, get all tags/taxadist from db"); 
            RepGenDataWriter repGenData=new RepGenSQLite(inputDB());

            Map <Tag, TaxaDistribution> tagTaxaMap = repGenData.getAllTagsTaxaMap();
            int tagcount = 0;

            int processedTags = 0;
            System.out.println("TIme to get all tags with taxa from db: " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
            time = System.nanoTime();
            Multimap<Tag,TagCorrelationInfo> tagTagCorrelations = Multimaps.synchronizedMultimap(HashMultimap.<Tag,TagCorrelationInfo>create());
            System.out.println("\nStart processing tag correlations.  Number of tags in db: " + tagTaxaMap.keySet().size());
            for (Tag tag1 : tagTaxaMap.keySet()) {
                tagcount++;
                processedTags++;
                // get dist for each taxa
                TaxaDistribution tag1TD = tagTaxaMap.get(tag1);
                if (tag1TD == null) {
                    ((RepGenSQLite)repGenData).close();
                    System.out.println("GetTagTaxaDist: got null tagTD at tagcount " + tagcount);
                    return null; // But this should return an error?
                }
                int[] depths1 = tag1TD.depths(); // gives us the depths for each taxon
                // This needs to be doubles for Pearson
                double[] ddepths1 = new double[depths1.length];
                
                // Apparently no shorter method of casting int to double 
                for (int idx = 0; idx < depths1.length; idx++) {
                    ddepths1[idx] = (double)depths1[idx];
                }

                double[] depthsPrime1 = new double[depths1.length];
                for (int idx = 0; idx < depthsPrime1.length; idx++) {
                    // boolean - presence/absence
                    if (ddepths1[idx] > 0) depthsPrime1[idx] = 1;
                    else depthsPrime1[idx] = 0;
                }
                tagTaxaMap.keySet().parallelStream().forEach(tag2 -> {
                    // don't skip if tag1= tag2?
                    TaxaDistribution tag2TD = tagTaxaMap.get(tag2);
                    if (tag2TD == null) {                       
                        System.out.println("GetTagTaxaDist: got null tagTD for sequence " + tag2.sequence());
                        return ;
                    }
                    // I need doubles below !!
                    int[] depths2 = tag2TD.depths(); // gives us the depths for each taxon
                    double[] ddepths2 = new double[depths2.length];
                    // Apparently no shorter method of casting int to double 
                    for (int idx = 0; idx < depths2.length; idx++) {
                        ddepths2[idx] = (double)depths2[idx];
                    }
                    double[] depthsPrime2 = new double[depths2.length];
                    for (int idx = 0; idx < depthsPrime1.length; idx++) {
                        // boolean - presence/absence
                        if (ddepths2[idx] > 0) depthsPrime2[idx] = 1;
                        else depthsPrime2[idx] = 0;
                    }
                    // From analysis.association.GenomicSelectionPlugin.java:
                    PearsonsCorrelation Pearsons = new PearsonsCorrelation();
                    double p1 = Pearsons.correlation(ddepths1,ddepths2);
                    
                    SpearmansCorrelation Spearmans = new SpearmansCorrelation();
                    double spearval = Spearmans.correlation(ddepths1, ddepths2);
                    
                    double p2 = Pearsons.correlation(depthsPrime1, depthsPrime2);
                    
                    // TODO LCJ Need to create r2 - This is from analysis.popgen.LinkageDisequilibrium
                    // Need to understand how to munge the data into the form it wants, how to 
                    // create genotype table.
                    
                    
                    TagCorrelationInfo tci = new TagCorrelationInfo(tag2,p1,spearval,p2,0);
                    tagTagCorrelations.put(tag1, tci);
                });
                
                if (tagcount > 1000) {
                    System.out.println("Finished processing " + processedTags + " tags");
                    tagcount = 0;
                }
            
            }
            if (tagcount > 0 ) {
                System.out.println("Finished processing last tags");
            }
            System.out.println("Total number of tags processed: " + processedTags);

            
            System.out.println("Size of tagTagCorrelations: " + tagTagCorrelations.size() + ", num tags in db: " 
              + tagTaxaMap.keySet().size() + ", num Tags in tagTagCorrelations map: " + tagTagCorrelations.keySet().size());
            
            // LCJ - this is debug
            Set<Tag> tags = tagTaxaMap.keySet();
            List<Tag> tagsAsList = new ArrayList<Tag>(tags);
            Collection<TagCorrelationInfo> tci = tagTagCorrelations.get(tagsAsList.get(0));
            System.out.println("LCJ - size of first tags correlation set: " + tci.size() + "\n");
            // end debug
            
            // Load to database     
            System.out.println("Loading correlation matrix to the database ...");
            repGenData.putTagTagCorrelationMatrix(tagTagCorrelations);
            ((RepGenSQLite)repGenData).close();
            
            System.out.println("Full RepGenLDAnalysis Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        } catch (Exception exc) {
            System.out.println("RepGenLDAnalysis:process_data:  processing error");
            exc.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }
    
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }
    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }
}
