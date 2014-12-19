/*
 * ProductionSNPCallerPlugin
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.atomic.LongAdder;

import javax.swing.ImageIcon;

import net.maizegenetics.analysis.gbs.Barcode;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.SimpleAllele;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.dna.snp.genotypecall.GenotypeMergeRule;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.dna.tag.TagData;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

/**
 * This plugin converts all of the fastq (and/or qseq) files in the input folder
 * and keyfile to genotypes and adds these to a genotype file in HDF5 format.
 *
 * We refer to this step as the "Production Pipeline".
 *
 * The output format is HDF5 genotypes with allelic depth stored. SNP calling is
 * quantitative with the option of using either the Glaubitz/Buckler binomial
 * method (pHet/pErr > 1 = het) (=default), or the Stacks method.
 *
 * Merging of samples with the same LibraryPrepID is handled by
 * GenotypeTableBuilder.addTaxon(), with the genotypes re-called based upon the
 * new depths. Therefore, if you want to keep adding genotypes to the same
 * target HDF5 file in subsequent runs, use the -ko (keep open) option so that
 * the output GenotypeTableBuilder will be mutable, using closeUnfinished()
 * rather than build().
 *
 * If the target output HDF5 GenotypeTable file doesn't exist, it will be
 * created.
 *
 * Each taxon in the HDF5 file is named "ShortName:LibraryPrepID" and is
 * annotated with "Flowcell_Lanes" (=source seq data for current genotype).
 *
 * Requires a TOPM with variants added from a previous "Discovery Pipeline" run.
 * In binary tagDataReader or HDF5 format (TOPMInterface).
 *
 * TODO add the Stacks likelihood method to BasicGenotypeMergeRule
 *
 * @author Ed Buckler
 * @author Jeff Glaubitz
 */
public class ProductionSNPCallerPluginV2 extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ProductionSNPCallerPluginV2.class);

    private PluginParameter<String> myInputDir = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Directory").required(true).inDir()
            .description("Input directory containing fastq AND/OR qseq files.").build();
    private PluginParameter<String> myKeyFile = new PluginParameter.Builder<>("k", null, String.class).guiName("Key File").required(true).inFile()
            .description("Key file listing barcodes distinguishing the samples").build();
    private PluginParameter<String> myEnzyme = new PluginParameter.Builder<>("e", null, String.class).guiName("Enzyme").required(true)
            .description("Enzyme used to create the GBS library").build();
    private PluginParameter<String> myInputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Input GBS Database").required(true).inFile()
            .description("Input Database file if using SQLite").build();
    private PluginParameter<String> myOutputGenotypes = new PluginParameter.Builder<>("o", null, String.class).guiName("Output HDF5 Genotypes File").required(true).outFile()
            .description("Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)").build();
    private PluginParameter<Double> myAveSeqErrorRate = new PluginParameter.Builder<>("eR", 0.01, Double.class).guiName("Ave Seq Error Rate")
            .description("Average sequencing error rate per base (used to decide between heterozygous and homozygous calls)").build();
    private PluginParameter<Integer> myMaxDivergence = new PluginParameter.Builder<>("d", 0, Integer.class).guiName("Max Divergence")
            .description("Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only)").build();
    private PluginParameter<Boolean> myKeepGenotypesOpen = new PluginParameter.Builder<>("ko", false, Boolean.class).guiName("Keep Genotypes Open")
            .description("Keep hdf5 genotypes open for future runs that add more taxa or more depth").build();
    private PluginParameter<Boolean> myNoDepthOutput = new PluginParameter.Builder<>("ndo", false, Boolean.class).guiName("No Depth to Output")
            .description("No depth output: do not write depths to the output hdf5 genotypes file").build();
    //private PluginParameter<Boolean> myStacksLikelihood = new PluginParameter.Builder<>("sL", false, Boolean.class).guiName("Use Stacks Likelihood")
    //        .description("Use STACKS likelihood method to call heterozygotes (default: use tasselGBS likelihood ratio method)").build();

    private String myOutputDir = null;
    private TagData tagDataReader = null;
    Multimap<Taxon,Tag> tagCntMap=Multimaps.synchronizedMultimap(ArrayListMultimap.create(384, 500_000));
    private Set<String> seqFilesInKeyAndDir = new TreeSet<>(); // fastq (or qseq) file names present in input directory that have a "Flowcell_Lane" in the key file

    private GenotypeTableBuilder genos = null; //output genotype table
    private PositionList myPositionList = null;

    //Documentation of read depth per sample (one recored per replicate)
    private Map<String, Integer> rawReadCountsForFullSampleName = new TreeMap<>();
    private Map<String, Integer> matchedReadCountsForFullSampleName = new TreeMap<>();

    private GenotypeMergeRule genoMergeRule = null;

    public ProductionSNPCallerPluginV2() {
        super(null, false);
    }

    public ProductionSNPCallerPluginV2(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public void postProcessParameters() {
        try {
            myOutputDir = (new File(outputHDF5GenotypesFile())).getCanonicalFile().getParent();
        } catch (IOException e) {
            throw new IllegalStateException("Problem resolving output directory:" + e);
        }
        genoMergeRule = new BasicGenotypeMergeRule(aveSeqErrorRate());
    }

    @Override
    public DataSet processData(DataSet input) {
        Path keyPath= Paths.get(keyFile()).toAbsolutePath();
        List<Path> inputSeqFiles= DirectoryCrawler.listPaths(GBSSeqToTagDBPlugin.inputFileGlob, Paths.get(myInputDir.value()).toAbsolutePath());
        if(inputSeqFiles.isEmpty()) {
            myLogger.warn("No files matching:"+GBSSeqToTagDBPlugin.inputFileGlob);
            return null;
        }
        tagDataReader =new TagDataSQLite(myInputDB.value());
        TaxaList masterTaxaList= TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), GBSSeqToTagDBPlugin.sampleNameField, new HashMap<>(), true);
        //todo perhaps subset the masterTaxaList based on the files in there, but it seems like it will all be figure out.
        Map<Tag,Tag> canonicalTag=new HashMap<>();  //canonicalize them OR eventually we will use a Trie
        tagDataReader.getTags().stream().forEach(t -> canonicalTag.put(t,t));
        inputSeqFiles.parallelStream()
                .forEach(inputSeqFile -> {
                    processFastQFile(masterTaxaList,keyPath, inputSeqFile, enzyme(),canonicalTag,64);
                });

        final PositionList positionList=tagDataReader.getSNPPositions();

        final Multimap<Tag,AlleleWithPosIndex> tagsToIndex=ArrayListMultimap.create();
        tagDataReader.getAlleleMap().entries().stream()
                .forEach(e -> {
                    int posIndex=positionList.indexOf(e.getValue().position());
                    tagsToIndex.put(e.getKey(),new AlleleWithPosIndex(e.getValue(),posIndex));
                });
        GenotypeTableBuilder gtb=setUpGenotypeTableBuilder(outputHDF5GenotypesFile(),positionList, genoMergeRule);
        tagCntMap.asMap().entrySet().stream()
                .map(e -> callGenotypes(e.getKey(), e.getValue(), tagsToIndex, positionList, genoMergeRule))
                .forEach(e -> gtb.addTaxon(e.x,e.y));
                //.forEach(e -> System.out.println(e.x.getName()+ Arrays.toString(Arrays.copyOfRange(e.y,0,10))));  //change this to a collector insert someplace


        if (keepGenotypesOpen()) {
            gtb.closeUnfinished();
        } else {
            gtb.build();
        }
        writeReadsPerSampleReports();
        return null;
    }

    private static Tuple<Taxon,byte[]> callGenotypes(Taxon taxon, Collection<Tag> tags, Multimap<Tag,AlleleWithPosIndex> tagsToIndex,
                                 PositionList positionList, GenotypeMergeRule genoMergeRule) {
        int[][] alleleDepths = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][positionList.numberOfSites()];
        tags.stream().map(t -> tagsToIndex.get(t)).flatMap(c -> c.stream())
                .forEach(a -> alleleDepths[a.allele()][a.positionIndex()]++);
        return new Tuple<>(taxon,resolveGenosForTaxon(alleleDepths, genoMergeRule));
    }

    private class AlleleWithPosIndex extends SimpleAllele {
        private int positionIndex;

        private AlleleWithPosIndex(Allele myAllele, int positionIndex) {
            super(myAllele.allele(), myAllele.position());
            this.positionIndex=positionIndex;
        }

        public int positionIndex() {
            return positionIndex;
        }
    }

    private class CountOfReadQuality {
        LongAdder allReads=new LongAdder();
        LongAdder goodBarcodedReads=new LongAdder();
        LongAdder goodMatched=new LongAdder();
        LongAdder perfectMatches=new LongAdder();
        LongAdder imperfectMatches=new LongAdder();
        LongAdder singleImperfectMatches=new LongAdder();
    }

    private void processFastQFile(TaxaList masterTaxaList, Path keyPath, Path fastQPath, String enzymeName,
                                  Map<Tag,Tag> canonicalTags, int preferredTagLength) {
        TaxaList tl=GBSSeqToTagDBPlugin.getLaneAnnotatedTaxaList(keyPath, fastQPath);
        BarcodeTrie barcodeTrie=GBSSeqToTagDBPlugin.initializeBarcodeTrie(tl, masterTaxaList, new GBSEnzyme(enzymeName));
        processFastQ(fastQPath,barcodeTrie,canonicalTags,preferredTagLength);
    }

    private void processFastQ(Path fastqFile, BarcodeTrie barcodeTrie, Map<Tag,Tag> canonicalTags, int preferredTagLength) {
        int allReads=0, goodBarcodedReads = 0;
        try {
            BufferedReader br = Utils.getBufferedReader(fastqFile.toString(), 1 << 22);
            long time=System.nanoTime();
            String[] seqAndQual;
            while ((seqAndQual=GBSSeqToTagDBPlugin.readFastQBlock(br, allReads)) != null) {
                allReads++;
                //After quality score is read, decode barcode using the current sequence & quality  score
                Barcode barcode=barcodeTrie.longestPrefix(seqAndQual[0]);
                if(barcode==null) continue;
                Tag tag= TagBuilder.instance(seqAndQual[0].substring(barcode.getBarLength(), barcode.getBarLength() + preferredTagLength)).build();
                if(tag==null) continue;   //null occurs when any base was not A, C, G, T
                goodBarcodedReads++;
                Tag canonicalTag=canonicalTags.get(tag);
                if(canonicalTag!=null) tagCntMap.put(barcode.getTaxon(),canonicalTag);
                if (allReads % 1000000 == 0) {
                    myLogger.info("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads
                            + " rate:" + (System.nanoTime()-time)/allReads +" ns/read");
                }
            }
            myLogger.info("Total number of reads in lane=" + allReads);
            myLogger.info("Total number of good barcoded reads=" + goodBarcodedReads);
            myLogger.info("Timing process (sorting, collapsing, and writing TagCount to file).");
            myLogger.info("Process took " + (System.nanoTime() - time)/1e6 + " milliseconds.");
            br.close();
        } catch (Exception e) {
            myLogger.error("Good Barcodes Read: " + goodBarcodedReads);
            e.printStackTrace();
        }
    }

    private void reportProgress(int[] counters, long readSeqReadTime, long ifRRNotNullTime) {
        myLogger.info(
                "totalReads:" + counters[0]
                + "  goodBarcodedReads:" + counters[1]
                + "  goodMatchedToTOPM:" + counters[2]
                //            + "  perfectMatches:" + counters[3]
                //            + "  nearMatches:" + counters[4]
                //            + "  uniqueNearMatches:" + counters[5]
                + "  cumulReadSequenceTime: " + ((double) (readSeqReadTime) / 1_000_000_000.0) + " sec"
                + "  cumulProcessSequenceTime: " + ((double) (ifRRNotNullTime) / 1_000_000_000.0) + " sec"
        );
    }

    private void reportTotals(Path fileName, int[] counters, int nFilesProcessed) {
        myLogger.info("Total number of reads in lane=" + counters[0]);
        myLogger.info("Total number of good, barcoded reads=" + counters[1]);
        myLogger.info("Total number of good, barcoded reads matched to the TOPM=" + counters[2]);
        myLogger.info("Finished reading " + nFilesProcessed + " of " + seqFilesInKeyAndDir.size() + " sequence files: " + fileName + "\n");
    }


    private static GenotypeTableBuilder setUpGenotypeTableBuilder(String aHDF5File, PositionList positionList, GenotypeMergeRule mergeRule) {
        File hdf5File = new File(aHDF5File);
        if (hdf5File.exists()) {
            myLogger.info("\nGenotypes will be added to existing HDF5 file:\n  " + aHDF5File + "\n");
            return GenotypeTableBuilder.mergeTaxaIncremental(aHDF5File, mergeRule);
        } else {
            myLogger.info("\nThe target HDF5 file:\n  " + aHDF5File
                    + "\ndoes not exist. A new HDF5 file of that name will be created \nto hold the genotypes from this run.");
            return GenotypeTableBuilder.getTaxaIncrementalWithMerging(aHDF5File, positionList, mergeRule);
        }
    }

    /**
     * Gets an ArrayList of taxa, each named "Sample:LibraryPrepID" and annotated
     * with Flowcell_Lane as well as all other annotations in the key file, for 
     * the corresponding fastq file.
     *
     * @return ArrayList<Taxon>
     */
//    private ArrayList<Taxon> getHDF5Taxa(int fileNum) {
//        String currFlowcellLane = seqFileNameToFlowcellLane.get(myRawSeqFileNames[fileNum]);
//        String[] flowcellLane = currFlowcellLane.split("_");
//        TaxaList annoTL = TaxaListIOUtils.readTaxaAnnotationFile(keyFile(), "LibraryPrepID",
//                ImmutableMap.of("Flowcell", flowcellLane[0], "Lane", flowcellLane[1]), false);
//        ArrayList<Taxon> taxaAL = new ArrayList();
//        for (Taxon tax : annoTL) {
//            String newName = tax.getTextAnnotation("Sample").length == 0 ? tax.getTextAnnotation("DNASample")[0] : tax.getTextAnnotation("Sample")[0];
//            String libPrepID = tax.getName();
//            newName += ":" + libPrepID;
//            Taxon gbsTaxon = new Taxon.Builder(tax).name(newName).addAnno("Flowcell_Lane", currFlowcellLane)
//                    .addAnno("LibraryPrepID", libPrepID).addAnno("Status", "private").build();
//            taxaAL.add(gbsTaxon);
//        }
//        return taxaAL;
//    }
//


//    private int findBestImperfectMatch(long[] read, int[] counters) {
//        // this method is not ready for prime time -- to resolve a tie, it currently chooses a random tag out of the tied tags
//        int tagIndex = -1;
//        TagMatchFinder tmf = new TagMatchFinder(tagDataReader);
//        TreeMap<Integer, Integer> bestHitsAndDiv = tmf.findMatchesWithIntLengthWords(read, maxDivergence(), true);
//        if (bestHitsAndDiv.size() > 0) {
//            counters[4]++; // imperfectMatches
//            if (bestHitsAndDiv.size() == 1) {
//                counters[5]++; // singleImperfectMatches
//            }
//            tagIndex = bestHitsAndDiv.firstKey();  // a random tag (firstKey) chosen to resolve the tie = suboptimal behavior
//        }
//        return tagIndex;
//    }
//
//    private void incrementDepthForTagVariants(int tagIndex, int[][] alleleDepths, int increment) {
//        int chromosome = tagDataReader.getChromosome(tagIndex);
//        if (chromosome == TOPMInterface.INT_MISSING) {
//            return;
//        }
//        int startPos = tagDataReader.getStartPosition(tagIndex);
//        for (int variant = 0; variant < tagDataReader.getMaxNumVariants(); variant++) {
//            byte newBase = tagDataReader.getVariantDef(tagIndex, variant);
//            if ((newBase == TOPMInterface.BYTE_MISSING) || (newBase == GenotypeTable.UNKNOWN_ALLELE)) {
//                continue;
//            }
//            int offset = tagDataReader.getVariantPosOff(tagIndex, variant);
//            int pos = startPos + offset;
////            int currSite = genos.getSiteOfPhysicalPosition(pos, locus);
//            int currSite = positionToSite.get(chromosome, pos);
//            if (currSite < 0) {
//                continue;
//            }
//            alleleDepths[newBase][currSite] += increment;
//        }
//    }

//    private void callGenotypes() {
//        myLogger.info("\nCalling genotypes...");
//        for (int currTaxonIndex = 0; currTaxonIndex < obsTagsForEachTaxon.length; currTaxonIndex++) {
//            IntArrayList currTagList = obsTagsForEachTaxon[currTaxonIndex];
//            currTagList.sort();
//            int[][] alleleDepths = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][myPositionList.numberOfSites()];
//            int prevTag = currTagList.getQuick(0);
//            int currInc = 0;
//            for (int t = 0; t < currTagList.size(); t++) {
//                int tag = currTagList.getQuick(t);
//                if (tag == prevTag) {
//                    currInc++;
//                } else {
//                    incrementDepthForTagVariants(prevTag, alleleDepths, currInc);
//                    prevTag = tag;
//                    currInc = 1;
//                }
//            }
//            incrementDepthForTagVariants(prevTag, alleleDepths, currInc);
//            byte[][] byteDepths = AlleleDepthUtil.depthIntToByte(alleleDepths);
//            byte[] taxonGenos = resolveGenosForTaxon(byteDepths);
//            if (noDepthToOutput()) {
//                genos.addTaxon(taxaList.get(currTaxonIndex), taxonGenos, null);
//            } else {
//                genos.addTaxon(taxaList.get(currTaxonIndex), taxonGenos, byteDepths);
//            }
//            myLogger.info("  finished calling genotypes for " + taxaList.get(currTaxonIndex).getName());
//        }
//        myLogger.info("Finished calling genotypes for " + obsTagsForEachTaxon.length + " taxa\n");
//    }

    private static byte[] resolveGenosForTaxon(int[][] depthsForTaxon, GenotypeMergeRule genoMergeRule) {
        int nAlleles = depthsForTaxon.length;
        int[] depthsAtSite = new int[nAlleles];
        int nSites = depthsForTaxon[0].length;
        byte[] genos = new byte[nSites];
        for (int site = 0; site < nSites; site++) {
            for (int allele = 0; allele < nAlleles; allele++) {
                depthsAtSite[allele] = depthsForTaxon[allele][site];
            }
            genos[site] = genoMergeRule.callBasedOnDepth(depthsAtSite);
        }
        return genos;
    }

    private void writeReadsPerSampleReports() {
        myLogger.info("\nWriting ReadsPerSample log file...");
        String outFileS = myOutputDir + File.separator + (new File(keyFile())).getName();
        outFileS = outFileS.replaceAll(".txt", "_ReadsPerSample.log");
        outFileS = outFileS.replaceAll("_key", "");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFileS))), 65536);
            bw.write("FullSampleName\tgoodBarcodedReads\tgoodReadsMatchedToTOPM\n");
            for (String fullSampleName : rawReadCountsForFullSampleName.keySet()) {
                bw.write(fullSampleName + "\t" + rawReadCountsForFullSampleName.get(fullSampleName) + "\t" + matchedReadCountsForFullSampleName.get(fullSampleName) + "\n");
            }
            bw.close();
        } catch (Exception e) {
            myLogger.error("Couldn't write to ReadsPerSample log file: " + e);
            e.printStackTrace();
            System.exit(1);
        }
        myLogger.info("   ...done\n");
    }

    private void printFileNameConventions(String actualFileName) {
        String message
                = "\n\n"
                + "Error in parsing file name:"
                + "\n   The raw sequence filename does not contain either 3, 4, or 5 underscore-delimited values."
                + "\n   Acceptable file naming conventions include the following (where FLOWCELL indicates the flowcell name and LANE is an integer):"
                + "\n       FLOWCELL_LANE_fastq.gz"
                + "\n       FLOWCELL_s_LANE_fastq.gz"
                + "\n       code_FLOWCELL_s_LANE_fastq.gz"
                + "\n       FLOWCELL_LANE_fastq.txt.gz"
                + "\n       FLOWCELL_s_LANE_fastq.txt.gz"
                + "\n       code_FLOWCELL_s_LANE_fastq.txt.gz"
                + "\n       FLOWCELL_LANE_qseq.txt.gz"
                + "\n       FLOWCELL_s_LANE_qseq.txt.gz"
                + "\n       code_FLOWCELL_s_LANE_qseq.txt.gz"
                + "\n"
                + "\n   Actual Filename: " + actualFileName
                + "\n\n";

        myLogger.error(message);
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Production SNP Caller";
    }

    @Override
    public String getToolTipText() {
        return "Production SNP Caller";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ProductionSNPCallerPluginV2.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TagData runPlugin(DataSet input) {
        return (TagData) performFunction(input).getData(0).getData();
    }

    /**
     * Input directory containing fastq AND/OR qseq files.
     *
     * @return Input Directory
     */
    public String inputDirectory() {
        return myInputDir.value();
    }

    /**
     * Set Input Directory. Input directory containing fastq
     * AND/OR qseq files.
     *
     * @param value Input Directory
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 inputDirectory(String value) {
        myInputDir = new PluginParameter<>(myInputDir, value);
        return this;
    }

    /**
     * Key file listing barcodes distinguishing the samples
     *
     * @return Key File
     */
    public String keyFile() {
        return myKeyFile.value();
    }

    /**
     * Set Key File. Key file listing barcodes distinguishing
     * the samples
     *
     * @param value Key File
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 keyFile(String value) {
        myKeyFile = new PluginParameter<>(myKeyFile, value);
        return this;
    }

    /**
     * Enzyme used to create the GBS library
     *
     * @return Enzyme
     */
    public String enzyme() {
        return myEnzyme.value();
    }

    /**
     * Set Enzyme. Enzyme used to create the GBS library
     *
     * @param value Enzyme
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 enzyme(String value) {
        myEnzyme = new PluginParameter<>(myEnzyme, value);
        return this;
    }

    /**
     * Input Database file if using SQLite
     *
     * @return Input GBS Database
     */
    public String inputGBSDatabase() {
        return myInputDB.value();
    }

    /**
     * Set Input GBS Database. Input Database file if using
     * SQLite
     *
     * @param value Input GBS Database
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 inputGBSDatabase(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
        return this;
    }

    /**
     * Output (target) HDF5 genotypes file to add new genotypes
     * to (new file created if it doesn't exist)
     *
     * @return Output HDF5 Genotypes File
     */
    public String outputHDF5GenotypesFile() {
        return myOutputGenotypes.value();
    }

    /**
     * Set Output HDF5 Genotypes File. Output (target) HDF5
     * genotypes file to add new genotypes to (new file created
     * if it doesn't exist)
     *
     * @param value Output HDF5 Genotypes File
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 outputHDF5GenotypesFile(String value) {
        myOutputGenotypes = new PluginParameter<>(myOutputGenotypes, value);
        return this;
    }

    /**
     * Average sequencing error rate per base (used to decide
     * between heterozygous and homozygous calls)
     *
     * @return Ave Seq Error Rate
     */
    public Double aveSeqErrorRate() {
        return myAveSeqErrorRate.value();
    }

    /**
     * Set Ave Seq Error Rate. Average sequencing error rate
     * per base (used to decide between heterozygous and homozygous
     * calls)
     *
     * @param value Ave Seq Error Rate
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 aveSeqErrorRate(Double value) {
        myAveSeqErrorRate = new PluginParameter<>(myAveSeqErrorRate, value);
        return this;
    }

    /**
     * Maximum divergence (edit distance) between new read
     * and previously mapped read (Default: 0 = perfect matches
     * only)
     *
     * @return Max Divergence
     */
    public Integer maxDivergence() {
        return myMaxDivergence.value();
    }

    /**
     * Set Max Divergence. Maximum divergence (edit distance)
     * between new read and previously mapped read (Default:
     * 0 = perfect matches only)
     *
     * @param value Max Divergence
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 maxDivergence(Integer value) {
        myMaxDivergence = new PluginParameter<>(myMaxDivergence, value);
        return this;
    }

    /**
     * Keep hdf5 genotypes open for future runs that add more
     * taxa or more depth
     *
     * @return Keep Genotypes Open
     */
    public Boolean keepGenotypesOpen() {
        return myKeepGenotypesOpen.value();
    }

    /**
     * Set Keep Genotypes Open. Keep hdf5 genotypes open for
     * future runs that add more taxa or more depth
     *
     * @param value Keep Genotypes Open
     *
     * @return this plugin
     */
    public ProductionSNPCallerPluginV2 keepGenotypesOpen(Boolean value) {
        myKeepGenotypesOpen = new PluginParameter<>(myKeepGenotypesOpen, value);
        return this;
    }

    /**
     * No depth output: do not write depths to the output
     * hdf5 genotypes file
     *
     * @return No Depth to Output
     */
    public Boolean noDepthToOutput() {
        return myNoDepthOutput.value();
    }


    /**
     * Use STACKS likelihood method to call heterozygotes (default: use
     * tasselGBS likelihood ratio method)
     *
     * @return Use Stacks Likelihood
     */
    //public Boolean useStacksLikelihood() {
    //    return myStacksLikelihood.value();
    //}
    /**
     * Set Use Stacks Likelihood. Use STACKS likelihood method to call
     * heterozygotes (default: use tasselGBS likelihood ratio method)
     *
     * @param value Use Stacks Likelihood
     *
     * @return this plugin
     */
    //public ProductionSNPCallerPlugin useStacksLikelihood(Boolean value) {
    //    myStacksLikelihood = new PluginParameter<>(myStacksLikelihood, value);
    //    return this;
    //}
}
