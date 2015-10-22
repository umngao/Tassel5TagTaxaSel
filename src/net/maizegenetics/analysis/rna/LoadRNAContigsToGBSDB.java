package net.maizegenetics.analysis.rna;



import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;


/**
 * Develops a discovery TBT file from a set of GBS sequence files.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 *
 * Originally the reference throughout was to "tag". This is being changed
 * to "kmer" as the pipeline is a kmer alignment process.
 *
 * @author Ed Buckler
 */
public class LoadRNAContigsToGBSDB extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(LoadRNAContigsToGBSDB.class);

    private PluginParameter<String> myContigFile = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Fasta Contig File").required(true).inFile()
            .description("Input file containing contigs in fasta format.\n").build();
    private PluginParameter<String> myOutputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Output Database File").required(true).outFile()
            .description("Output Database File").build();
    private PluginParameter<Boolean> myDeleteOldData = new PluginParameter.Builder<Boolean>("deleteOldData",false,Boolean.class).guiName("Delete Old Data")
            .description("Delete existing SNP quality data from db tables").build();

    public LoadRNAContigsToGBSDB() {
        super(null, false);
    }

    public LoadRNAContigsToGBSDB(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        super.preProcessParameters(input);
        if (Files.exists(Paths.get(outputDB()))) {
            if (deleteOldData()) {
                try {
                    Files.delete(Paths.get(myOutputDB.value()));
                } catch (Exception exc) {
                    System.out.println("Error when trying to delete database file: " + myOutputDB.value());
                    System.out.println("File delete error: " + exc.getMessage());
                }
            }
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        long lineNumber=0;

        try {
            TagDataWriter tdw =new TagDataSQLite(myOutputDB.value());
            BufferedReader br = Utils.getBufferedReader(contigFile(), 1 << 22);
            String line;
            String header=null;
            StringBuilder seq=new StringBuilder();
            Map<Tag,String> contigNameMap=new HashMap<>(100_000);
            while ((line = br.readLine())!= null) {
                line = line.trim();
                lineNumber++;
                if (line.startsWith(">")) {
                    if(header!=null) {
                        contigNameMap.put(TagBuilder.instance(seq.toString()).build(), header);
                    }
                    //reset to new sequence
                    header=line;
                    seq=new StringBuilder();
                } else {
                    seq.append(line);
                }
            }
            if(header!=null) {
                contigNameMap.put(TagBuilder.instance(seq.toString()).build(),header);
            }
            tdw.putAllNamesTag(contigNameMap);  //add map to databse
            ((TagDataSQLite)tdw).close();

        } catch(Exception ioe) {
            System.err.println("Error in line number "+lineNumber);
            ioe.printStackTrace();
        }
        return new DataSet(new Datum("OutputDatabase",outputDB(),""),this);
    }




    private void processFastQFile(Path fastQPath) throws StringIndexOutOfBoundsException {
//        ArrayList<Taxon> tl=GBSUtils.getLaneAnnotatedTaxaList(keyPath, fastQPath);
//        if (tl.size() == 0) return;
//        BarcodeTrie barcodeTrie=GBSUtils.initializeBarcodeTrie(tl, masterTaxaList, new GBSEnzyme(enzymeName));
//        try {
//            processFastQ(fastQPath,barcodeTrie,masterTaxaList,masterTagTaxaMap,preferredTagLength,minQuality);
//        } catch (StringIndexOutOfBoundsException oobe) {
//            throw oobe; // Let processData() handle it - we want to stop processing on this error
//        }
    }




// The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
//     public static void main(String[] args) {
//         GeneratePluginCode.generate(LoadRNAContigsToGBSDB.class);
//     }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public String runPlugin(DataSet input) {
        return (String) performFunction(input).getData(0).getData();
    }

    /**
     * Input file containing contigs in fasta format.
     *
     *
     * @return Input Fasta Contig File
     */
    public String contigFile() {
        return myContigFile.value();
    }

    /**
     * Set Input Fasta Contig File. Input file containing
     * contigs in fasta format.
     *
     *
     * @param value Input Fasta Contig File
     *
     * @return this plugin
     */
    public LoadRNAContigsToGBSDB contigFile(String value) {
        myContigFile = new PluginParameter<>(myContigFile, value);
        return this;
    }

    /**
     * Output Database File
     *
     * @return Output Database File
     */
    public String outputDB() {
        return myOutputDB.value();
    }

    /**
     * Set Output Database File. Output Database File
     *
     * @param value Output Database File
     *
     * @return this plugin
     */
    public LoadRNAContigsToGBSDB outputDB(String value) {
        myOutputDB = new PluginParameter<>(myOutputDB, value);
        return this;
    }

    /**
     * Delete existing SNP quality data from db tables
     *
     * @return Delete Old Data
     */
    public Boolean deleteOldData() {
        return myDeleteOldData.value();
    }

    /**
     * Set Delete Old Data. Delete existing SNP quality data
     * from db tables
     *
     * @param value Delete Old Data
     *
     * @return this plugin
     */
    public LoadRNAContigsToGBSDB deleteOldData(Boolean value) {
        myDeleteOldData = new PluginParameter<>(myDeleteOldData, value);
        return this;
    }


    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }

    @Override
    public String getToolTipText() {
        return null;
    }
}


