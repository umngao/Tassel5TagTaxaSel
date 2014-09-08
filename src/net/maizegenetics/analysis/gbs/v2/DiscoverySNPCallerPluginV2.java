/*
 * DiscoverySNPCallerPlugin
 */
package net.maizegenetics.analysis.gbs.v2;

import com.google.common.collect.*;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.util.ConcurrencyTools;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.IntFunction;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.maizegenetics.plugindef.PluginParameter;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.*;

/**
 * This class aligns tags at the same physical location against one another,
 * calls SNPs, and then outputs the SNPs to a HapMap file.
 *
 * It is multi-threaded, as there are substantial speed increases with it.
 *
 * @author Jeff Glaubitz
 * @author Ed Buckler
 */
public class DiscoverySNPCallerPluginV2 extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(DiscoverySNPCallerPluginV2.class);

    private PluginParameter<String> myInputTagsByTaxa = new PluginParameter.Builder<>("i", null, String.class).guiName("Input Tags by Taxa File").required(true).inFile()
            .description("Input TagsByTaxa file (if hdf5 format, use .hdf or .h5 extension)").build();
    private PluginParameter<String> myLogFile = new PluginParameter.Builder<>("log", null, String.class).guiName("Log File").outFile()
            .description("TagLocus log file name. (Default: TagLocusLog.txt)").build();
    private PluginParameter<Double> myMinMinorAlleleFreq = new PluginParameter.Builder<>("mnMAF", 0.01, Double.class).guiName("Min Minor Allele Freq")
            .description("Minimum minor allele frequency").build();
    private PluginParameter<Integer> myMinMinorAlleleCount = new PluginParameter.Builder<>("mnMAC", 10, Integer.class).guiName("Min Minor Allele Count")
            .description("Minimum minor allele count").build();
    private PluginParameter<Double> myMinLocusCoverage = new PluginParameter.Builder<>("mnLCov", 0.1, Double.class).guiName("Min Locus Coverage")
            .description("Minimum locus coverage (proportion of Taxa with a genotype)").build();
    private PluginParameter<String> myRefGenome = new PluginParameter.Builder<>("ref", null, String.class).guiName("Reference Genome File").inFile()
            .description("Path to reference genome in fasta format. Ensures that a tag from the reference genome is always included "
                    + "when the tags at a locus are aligned against each other to call SNPs. The reference allele for each site "
                    + "is then provided in the output HapMap files, under the taxon name \"REFERENCE_GENOME\" (first taxon). "
                    + "DEFAULT: Don't use reference genome.").build();
    private PluginParameter<Integer> myStartChr = new PluginParameter.Builder<>("sC", null, Integer.class).guiName("Start Chromosome").required(true)
            .description("Start Chromosome").build();
    private PluginParameter<Integer> myEndChr = new PluginParameter.Builder<>("eC", null, Integer.class).guiName("End Chromosome").required(true)
            .description("End Chromosome").build();
    private PluginParameter<Boolean> myIncludeRareAlleles = new PluginParameter.Builder<>("inclRare", false, Boolean.class).guiName("Include Rare Alleles")
            .description("Include the rare alleles at site (3 or 4th states)").build();
    private PluginParameter<Boolean> myIncludeGaps = new PluginParameter.Builder<>("inclGaps", false, Boolean.class).guiName("Include Gaps")
            .description("Include sites where major or minor allele is a GAP").build();
    private PluginParameter<Boolean> myCallBiSNPsWGap = new PluginParameter.Builder<>("callBiSNPsWGap", false, Boolean.class).guiName("Call Biallelic SNPs with Gap")
            .description("Include sites where the third allele is a GAP (mutually exclusive with inclGaps)").build();

   // private TagsOnPhysicalMap theTOPM = null;
    private int minTaxaWithLocus;
    private TagDataWriter theTBT = null;
    private boolean includeReference = false;
    private long[] refGenomeChr = null;
    private boolean customSNPLogging = true;  // a custom SNP log that collects useful info for filtering SNPs through machine learning criteria
//    private CustomSNPLog myCustomSNPLog = null;
    private boolean customFiltering = false;

    public DiscoverySNPCallerPluginV2() {
        super(null, false);
    }

    public DiscoverySNPCallerPluginV2(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        myLogger.info("Finding SNPs in " + inputTagsByTaxaFile() + ".");
        myLogger.info(String.format("StartChr:%d EndChr:%d %n", startChromosome(), endChromosome()));
        DataOutputStream locusLogDOS = openLocusLog(logFile());
        if (customSNPLogging) {
//            myCustomSNPLog = new CustomSNPLog(logFile());
        }
        for (int chr = startChromosome(); chr <= endChromosome(); chr++) {
            myLogger.info("\n\nProcessing chromosome " + chr + "...");
            if (includeReference) {
                refGenomeChr = readReferenceGenomeChr(referenceGenomeFile(), chr);
                if (refGenomeChr == null) {
                    myLogger.info("  WARNING: chromosome " + chr + " not found in the reference genome file. Skipping this chromosome.");
                    continue;
                }
            }
            discoverSNPsOnChromosome(chr, locusLogDOS);
            myLogger.info("Finished processing chromosome " + chr + "\n\n");
        }

        try {
            locusLogDOS.close();
        } catch (Exception e) {
            catchLocusLogException(e);
        }
        if (customSNPLogging) {
//            myCustomSNPLog.close();
        }
        ConcurrencyTools.shutdown();
        return null;
    }

    @Override
    public void postProcessParameters() {

        if (myInputTagsByTaxa.isEmpty()) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Input Tags by Taxa File not Set.");
        } else {
            theTBT=new TagDataSQLite(inputTagsByTaxaFile());
        }
        if (theTBT == null) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Problem reading Tags by Taxa File: " + inputTagsByTaxaFile());
        }
//
//        boolean loadBinary = (inputTOPMFile().endsWith(".txt")) ? false : true;
//        theTOPM = new TagsOnPhysicalMap(inputTOPMFile(), loadBinary);
//
//        if (myLogFile.isEmpty()) {
//            try {
//                File outDir = (new File(outputTOPMFile())).getCanonicalFile().getParentFile();
//                logFile(outDir.getCanonicalPath() + File.separator + "TagLocusLog.txt");
//            } catch (IOException e) {
//                throw new IllegalArgumentException("Problem creating the tagLocusLog file. Program aborted: "+e);
//            }
//        }
//
//        minTaxaWithLocus = (int) Math.round(theTBT.getTaxaCount() * minLocusCoverage());
//
        if (!myRefGenome.isEmpty()) {
            includeReference = true;
        }

        if (callBiallelicSNPsWithGap() && includeGaps()) {
            throw new IllegalArgumentException("The callBiSNPsWGap option is mutually exclusive with the inclGaps option.");
        }

        if (endChromosome() - startChromosome() < 0) {
            throw new IllegalArgumentException("The start chromosome is larger than the end chromosome.");
        }

        myLogger.info(String.format("minTaxaWithLocus:%d MinMAF:%g MinMAC:%d %n", minTaxaWithLocus, minMinorAlleleFreq(), minMinorAlleleCount()));
        myLogger.info(String.format("includeRare:%s includeGaps:%s %n", includeRareAlleles(), includeGaps()));
    }

    //
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(DiscoverySNPCallerPlugin.class);
    // }
    //
    /**
     * Input TagsByTaxa file (if hdf5 format, use .hdf or .h5 extension)
     *
     * @return Input Tags by Taxa File
     */
    public String inputTagsByTaxaFile() {
        return myInputTagsByTaxa.value();
    }

    /**
     * Set Input Tags by Taxa File. Input TagsByTaxa file (if hdf5 format, use
     * .hdf or .h5 extension)
     *
     * @param value Input Tags by Taxa File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 inputTagsByTaxaFile(String value) {
        myInputTagsByTaxa = new PluginParameter<>(myInputTagsByTaxa, value);
        return this;
    }

    /**
     * TagLocus log file name
     *
     * @return Log File
     */
    public String logFile() {
        return myLogFile.value();
    }

    /**
     * Set Log File. TagLocus log file name
     *
     * @param value Log File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 logFile(String value) {
        myLogFile = new PluginParameter<>(myLogFile, value);
        return this;
    }

    /**
     * Minimum minor allele frequency
     *
     * @return Min Minor Allele Freq
     */
    public Double minMinorAlleleFreq() {
        return myMinMinorAlleleFreq.value();
    }

    /**
     * Set Min Minor Allele Freq. Minimum minor allele frequency
     *
     * @param value Min Minor Allele Freq
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 minMinorAlleleFreq(Double value) {
        myMinMinorAlleleFreq = new PluginParameter<>(myMinMinorAlleleFreq, value);
        return this;
    }

    /**
     * Minimum minor allele count
     *
     * @return Min Minor Allele Count
     */
    public Integer minMinorAlleleCount() {
        return myMinMinorAlleleCount.value();
    }

    /**
     * Set Min Minor Allele Count. Minimum minor allele count
     *
     * @param value Min Minor Allele Count
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 minMinorAlleleCount(Integer value) {
        myMinMinorAlleleCount = new PluginParameter<>(myMinMinorAlleleCount, value);
        return this;
    }

    /**
     * Minimum locus coverage (proportion of Taxa with a genotype)
     *
     * @return Min Locus Coverage
     */
    public Double minLocusCoverage() {
        return myMinLocusCoverage.value();
    }

    /**
     * Set Min Locus Coverage. Minimum locus coverage (proportion of Taxa with a
     * genotype)
     *
     * @param value Min Locus Coverage
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 minLocusCoverage(Double value) {
        myMinLocusCoverage = new PluginParameter<>(myMinLocusCoverage, value);
        return this;
    }

    /**
     * Path to reference genome in fasta format. Ensures that a tag from the
     * reference genome is always included when the tags at a locus are aligned
     * against each other to call SNPs. The reference allele for each site is
     * then provided in the output HapMap files, under the taxon name
     * "REFERENCE_GENOME" (first taxon). DEFAULT: Don't use reference genome.
     *
     * @return Reference Genome File
     */
    public String referenceGenomeFile() {
        return myRefGenome.value();
    }

    /**
     * Set Reference Genome File. Path to reference genome in fasta format.
     * Ensures that a tag from the reference genome is always included when the
     * tags at a locus are aligned against each other to call SNPs. The
     * reference allele for each site is then provided in the output HapMap
     * files, under the taxon name "REFERENCE_GENOME" (first taxon). DEFAULT:
     * Don't use reference genome.
     *
     * @param value Reference Genome File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 referenceGenomeFile(String value) {
        myRefGenome = new PluginParameter<>(myRefGenome, value);
        return this;
    }

    /**
     * Include the rare alleles at site (3 or 4th states)
     *
     * @return Include Rare Alleles
     */
    public Boolean includeRareAlleles() {
        return myIncludeRareAlleles.value();
    }

    /**
     * Set Include Rare Alleles. Include the rare alleles at site (3 or 4th
     * states)
     *
     * @param value Include Rare Alleles
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 includeRareAlleles(Boolean value) {
        myIncludeRareAlleles = new PluginParameter<>(myIncludeRareAlleles, value);
        return this;
    }

    /**
     * Include sites where major or minor allele is a GAP
     *
     * @return Include Gaps
     */
    public Boolean includeGaps() {
        return myIncludeGaps.value();
    }

    /**
     * Set Include Gaps. Include sites where major or minor allele is a GAP
     *
     * @param value Include Gaps
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 includeGaps(Boolean value) {
        myIncludeGaps = new PluginParameter<>(myIncludeGaps, value);
        return this;
    }

    /**
     * Include sites where the third allele is a GAP (mutually exclusive with
     * inclGaps)
     *
     * @return Call Biallelic SNPs with Gap
     */
    public Boolean callBiallelicSNPsWithGap() {
        return myCallBiSNPsWGap.value();
    }

    /**
     * Set Call Biallelic SNPs with Gap. Include sites where the third allele is
     * a GAP (mutually exclusive with inclGaps)
     *
     * @param value Call Biallelic SNPs with Gap
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 callBiallelicSNPsWithGap(Boolean value) {
        myCallBiSNPsWGap = new PluginParameter<>(myCallBiSNPsWGap, value);
        return this;
    }

    /**
     * Start Chromosome
     *
     * @return Start Chromosome
     */
    public Integer startChromosome() {
        return myStartChr.value();
    }

    /**
     * Set Start Chromosome. Start Chromosome
     *
     * @param value Start Chromosome
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 startChromosome(Integer value) {
        myStartChr = new PluginParameter<>(myStartChr, value);
        return this;
    }

    /**
     * End Chromosome
     *
     * @return End Chromosome
     */
    public Integer endChromosome() {
        return myEndChr.value();
    }

    /**
     * Set End Chromosome. End Chromosome
     *
     * @param value End Chromosome
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 endChromosome(Integer value) {
        myEndChr = new PluginParameter<>(myEndChr, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Discovery SNP Caller";
    }

    @Override
    public String getToolTipText() {
        return "Discovery SNP Caller";
    }

    public void discoverSNPsOnChromosome(int targetChromo, DataOutputStream locusLogDOS) {
        int countLoci = 0;
        int siteCnt = 0;
        theTBT.getCutPositionTagTaxaMap(new Chromosome(""+targetChromo),-1,-1).entrySet().stream()
                .forEach(emp -> findAlleleByAlignment(emp.getKey(),emp.getValue()));
        myLogger.info("Number of marker sites recorded for chr" + targetChromo + ": " + siteCnt);
    }

    private void findAlleleByAlignment(Position position, Map<Tag,TaxaDistribution> tagTaxaMap) {
        if(tagTaxaMap.size()<2) {//homozygous
            //logRejectedTagLocus(theTAL, locusLogDOS);
            return;  //consider reporting homozygous
        }
        int taxaCoverage=tagTaxaMap.values().stream().mapToInt(TaxaDistribution::totalDepth).sum();  //todo this could be changed to taxa with tag
        if(taxaCoverage < myMinLocusCoverage.value()) {
            return;  //consider reporting low coverage
        }
        Map<Tag,String> alignedTags=alignTags(tagTaxaMap.keySet());
        Table<Position, Byte, List<TagTaxaDistribution>> tAlign=convertAlignmentToTagTable(alignedTags, tagTaxaMap,  position);
        //Filter tAlign by reference present
        //Filter tAlign by MAF & coverage
        //update the db
        //log information
    }

    static Map<Tag,String> alignTags(Collection<Tag> tags) {
        List<DNASequence> lst = tags.stream()
                .map(t -> {
                    DNASequence ds=new DNASequence(t.sequence());
                    ds.setUserCollection(ImmutableList.of(t));
                    return ds;
                })
                .collect(Collectors.toList());
        Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        ImmutableMap.Builder<Tag,String> result=new ImmutableMap.Builder<>();
        for (AlignedSequence<DNASequence, NucleotideCompound> compounds : profile) {
            ImmutableList tagList=(ImmutableList)compounds.getOriginalSequence().getUserCollection();
            result.put((Tag)tagList.get(0),compounds.getSequenceAsString());
        }

        return result.build();
    }

    /**
     * Converts at TagAlignment
     * @param alignedTags
     * @param tagTaxaDistMap
     * @param refStartPosition
     * @return
     */
    public static Table<Position, Byte, List<TagTaxaDistribution>> convertAlignmentToTagTable(Map<Tag,String> alignedTags,
                        Map<Tag,TaxaDistribution> tagTaxaDistMap, Position refStartPosition) {
        Table<Position, Byte, List<TagTaxaDistribution>> alignT= TreeBasedTable.create(); //These could be sorted by depth
        alignedTags.forEach((t,s) -> {
            TagTaxaDistribution td=new TagTaxaDistribution(t,tagTaxaDistMap.get(t));
            int start2=refStartPosition.getPosition();
            for (int i = 0; i < s.length(); i++) {
                byte allele=getNucleotideAlleleByte(s.charAt(i));
                Position p=new GeneralPosition.Builder(refStartPosition.getChromosome(),start2+i).build();
                List tdList=alignT.get(p,allele);
                if(tdList!=null) {tdList.add(td);}
                else {
                    List<TagTaxaDistribution> ttdL=new ArrayList<>();
                    ttdL.add(td);
                    alignT.put(p,allele, ttdL);
                }
            }
        });
        return alignT;
    }

    public static int tagDepthSum(List<TagTaxaDistribution> tagTaxaDistributions) {
        return tagTaxaDistributions.stream().mapToInt(t -> t.taxaDist().totalDepth()).sum();
    }

    public static int tagDepthSum(Map<Byte,List<TagTaxaDistribution>> tagTaxaDistributions) {
        return tagTaxaDistributions.values().stream().mapToInt(t -> tagDepthSum(t)).sum();
    }

    @Deprecated
    public static int[][] baseDepths(Map<Tag,String> alignedTags, Map<Tag,TaxaDistribution> tagTaxaDistMap) {
        int alignLength=alignedTags.values().stream().findFirst().map(String::length).get();
        int[][] alleleDepth=new int[alignLength][6];
        alignedTags.forEach((t,s) -> {
            TaxaDistribution td=tagTaxaDistMap.get(t);
            for (int i = 0; i < s.length(); i++) {
                int allele=NucleotideAlignmentConstants.getNucleotideAlleleByte(s.charAt(i));
                alleleDepth[i][allele]+=td.totalDepth();
            }
        });
        return alleleDepth;
    }

    public static List<Optional<Integer>> referencePositions(Position refStartPosition, Map<Tag,String> alignedTags){
        Tag refTag=alignedTags.keySet().stream()
                .filter(t -> t.isReference()).findFirst().orElseThrow(() -> new IllegalStateException("Reference not found"));
        AtomicInteger start=new AtomicInteger(refStartPosition.getPosition());
        List<Optional<Integer>> refPositionInts2=alignedTags.get(refTag).chars()
                .mapToObj(c -> (c == GAP_ALLELE_CHAR)?Optional.<Integer>empty():Optional.of(start.getAndIncrement()))
                .collect(Collectors.toList());
        return refPositionInts2;
    }

//    private void updateTOPM(TagRefAlignment myTAL, int variableSite, int position, int strand, byte[] alleles) {
//        for (int tg = 0; tg < myTAL.getSize(); tg++) {
//            int topmTagIndex = myTAL.getTOPMIndexOfTag(tg);
//            if (topmTagIndex == Integer.MIN_VALUE) {
//                continue; // skip the reference genome tag (which may not be in the TOPM)
//            }
//            byte baseToAdd = myTAL.getCallAtVariableSiteForTag(variableSite, tg);
//            boolean matched = false;
//            for (byte cb : alleles) {
//                if (baseToAdd == cb) {
//                    matched = true;
//                    break;
//                }
//            }
//            // so that all tags in the tagAlignment have the same corresponding variants in the TOPM, add a variant no matter what (set to missing if needed)
//            byte offset = (byte) (position - myTAL.getMinStartPosition());
//            if (!matched) {
//                baseToAdd = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
//            }
//            if (strand == -1) {
//                baseToAdd = NucleotideAlignmentConstants.getNucleotideComplement(baseToAdd);  // record everything relative to the plus strand
//            }
//            // convert from allele from 0-15 style to IUPAC ASCII character value (e.g., (byte) 'A') (maintains compatibility with Tassel3 TOPM)
//            baseToAdd = getIUPACAllele(baseToAdd);
//            theTOPM.addVariant(topmTagIndex, offset, baseToAdd);
//        }
//    }


    private DataOutputStream openLocusLog(String logFileName) {
        try {
            DataOutputStream locusLogDOS
                    = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(logFileName)), 65536));
            locusLogDOS.writeBytes(
                    "chr\tstart\tend\tstrand\ttotalbp\tnTags\tnReads\tnTaxaCovered\tminTaxaCovered\tstatus\tnVariableSites\tposVariableSites\tnVarSitesKept\tposVarSitesKept\trefTag?\tmaxTagLen\tminTagLen\n");
            return locusLogDOS;
        } catch (Exception e) {
            catchLocusLogException(e);
        }
        return null;
    }

//    private void logRejectedTagLocus(TagRefAlignment currTAL, DataOutputStream locusLogDOS) {
//        int start, end;
//        if (currTAL.getStrand() == -1) {
//            end = currTAL.getMinStartPosition();
//            start = currTAL.getMinStartPosition() - currTAL.getMaxTagLength() + 1;
//        } else {
//            start = currTAL.getMinStartPosition();
//            end = currTAL.getMinStartPosition() + currTAL.getMaxTagLength() - 1;
//        }
//        int totalbp = end - start + 1;
//        String status, refTag;
//        if (currTAL.getSize() == 1) {
//            status = "invariant\t0";
//            refTag = currTAL.getDivergenceOfTag(0) == 0 ? "1" : "0";
//        } else {
//            status = "tooFewTaxa\tNA";
//            boolean refTagFound = false;
//            int t = -1;
//            while (!refTagFound && t < currTAL.getSize() - 1) {
//                t++;
//                if (currTAL.getDivergenceOfTag(t) == 0) {
//                    refTagFound = true;
//                }
//            }
//            refTag = refTagFound ? "1" : "0";
//        }
//        try {
//            locusLogDOS.writeBytes(
//                    currTAL.getChromosome() + "\t"
//                    + start + "\t"
//                    + end + "\t"
//                    + currTAL.getStrand() + "\t"
//                    + totalbp + "\t"
//                    + currTAL.getSize() + "\t"
//                    + currTAL.getTotalNReads() + "\t"
//                    + currTAL.getNumberTaxaCovered() + "\t"
//                    + minTaxaWithLocus + "\t"
//                    + status + "\t"
//                    + "NA" + "\t"
//                    + "0" + "\t"
//                    + "NA" + "\t"
//                    + refTag + "\t"
//                    + currTAL.getMaxTagLength() + "\t"
//                    + currTAL.getMinTagLength() + "\n"
//            );
//        } catch (Exception e) {
//            catchLocusLogException(e);
//        }
//    }

    private void logAcceptedTagLocus(String locusLogRecord, DataOutputStream locusLogDOS) {
        try {
            locusLogDOS.writeBytes(locusLogRecord);
        } catch (Exception e) {
            catchLocusLogException(e);
        }
    }

    private void catchLocusLogException(Exception e) {
        System.out.println("ERROR: Unable to write to locus log file: " + e);
        e.printStackTrace();
        System.exit(1);
    }


    //TODO: this should be a sequence class.  It should also be rewritten with Java 8 & nio tools
    private long[] readReferenceGenomeChr(String inFileStr, int targetChr) {
        int nBases = getLengthOfReferenceGenomeChr(inFileStr, targetChr);
        if (nBases == 0) {
            return null;
        }
        int basesPerLong = BaseEncoder.chunkSize;
        int nLongs = (nBases % basesPerLong == 0) ? nBases / basesPerLong : (nBases / basesPerLong) + 1;
        long[] refGenomeChrAsLongs = new long[nLongs];
        myLogger.info("\n\nReading in the target chromosome " + targetChr + " from the reference genome fasta file: " + inFileStr);
        String temp = "Nothing has been read yet from the reference genome fasta file";
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(inFileStr)));
            StringBuilder currStrB = new StringBuilder();
            int currChr = Integer.MIN_VALUE, chunk = 0;
            while (br.ready()) {
                temp = br.readLine().trim();
                if (temp.startsWith(">")) {
                    if (chunk > 0) {
                        break;  // finished reading the targetChr (no need to read the rest of the file)
                    }
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    currChr = Integer.parseInt(chrS);  // don't need to catch exception because getLengthOfReferenceGenomeChr() would have caught it already
                    myLogger.info("Currently reading chromosome " + currChr + " (target chromosome = " + targetChr + ")");
                } else if (currChr == targetChr) {
                    currStrB.append(temp.replace("N", "A")); // BaseEncoder encodes sequences with N's as (long) -1
                    while (currStrB.length() >= basesPerLong) {
                        refGenomeChrAsLongs[chunk] = BaseEncoder.getLongFromSeq(currStrB.substring(0, basesPerLong));
                        currStrB = (currStrB.length() > basesPerLong) ? new StringBuilder(currStrB.substring(basesPerLong)) : new StringBuilder();
                        chunk++;
                        if (chunk % 1000000 == 0) {
                            myLogger.info(chunk + " chunks of " + basesPerLong + " bases read from the reference genome fasta file for chromosome " + targetChr);
                        }
                    }
                }
            }
            if (currStrB.length() > 0) {
                refGenomeChrAsLongs[chunk] = BaseEncoder.getLongFromSeq(currStrB.toString());
                chunk++;
            }
            myLogger.info("\n\nFinished reading target chromosome " + targetChr + " into a total of " + chunk + " " + basesPerLong + "bp chunks\n\n");
            if (chunk != nLongs) {
                throw new Exception("The number of 32 base chunks read (" + chunk + ") was not equal to the expected number (" + nLongs + ")");
            }
            br.close();
        } catch (Exception e) {
            myLogger.error("Exception caught while reading the reference genome fasta file at line.  Error=" + e);
            e.printStackTrace();
            System.exit(1);
        }
        return refGenomeChrAsLongs;
    }

    //TODO this should be a sequence class.  It should also be rewritten with Java 8 & nio tools
    private int getLengthOfReferenceGenomeChr(String inFileStr, int targetChr) {
        myLogger.info("\n\nDetermining the length (in bases) of target chromosome " + targetChr + " in the reference genome fasta file: " + inFileStr);
        String temp = "Nothing has been read yet from the reference genome fasta file";
        int line = 0, nBases = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(inFileStr)));
            int currChr = Integer.MIN_VALUE;
            while (br.ready()) {
                temp = br.readLine().trim();
                line++;
                if (line % 1000000 == 0) {
                    myLogger.info(line + " lines read from the reference genome fasta file");
                }
                if (temp.startsWith(">")) {
                    if (nBases > 0) {
                        break;  // finished reading the targetChr (no need to read the rest of the file)
                    }
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    try {
                        currChr = Integer.parseInt(chrS);
                    } catch (NumberFormatException e) {
                        myLogger.error("\n\nTagsToSNPByAlignment detected a non-numeric chromosome name in the reference genome sequence fasta file: " + chrS
                                + "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
                                + "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n");
                        System.exit(1);
                    }
                    myLogger.info("Currently reading chromosome " + currChr + " (target chromosome = " + targetChr + ")");
                } else if (currChr == targetChr) {
                    nBases += temp.length();
                }
            }
            if (nBases == 0) {
                throw new Exception("Target chromosome (" + targetChr + ") not found");
            }
            myLogger.info("The target chromosome " + targetChr + " is " + nBases + " bases long");
            br.close();
        } catch (Exception e) {
            if (nBases == 0) {
                myLogger.warn("Exception caught while reading the reference genome fasta file at line " + line + "\n   e=" + e + "\n   Skipping this chromosome...");
            } else {
                myLogger.error("Exception caught while reading the reference genome fasta file at line " + line + "\n   e=" + e);
                e.printStackTrace();
                System.exit(1);
            }
        }
        return nBases;
    }


    //TODO this should be a sequence class.  It should also be rewritten with Java 8 & nio tools
//    private String getRefSeqInRegion(TagRefAlignment theTAL) {
//        int basesPerLong = BaseEncoder.chunkSize;
//        int refSeqStartPosition = theTAL.getMinStartPosition() - 128;
//        int startIndex = Math.max((refSeqStartPosition / basesPerLong) - 1, 0);
//        int refSeqEndPosition = theTAL.getMaxStartPosition() + 128;
//        int endIndex = Math.min((refSeqEndPosition / basesPerLong) + 1, refGenomeChr.length - 1);
//        StringBuilder sb = new StringBuilder();
//        for (int i = startIndex; i <= endIndex; ++i) {
//            sb.append(BaseEncoder.getSequenceFromLong(refGenomeChr[i]));
//        }
//        theTAL.setMinStartPosition(startIndex * basesPerLong + 1);
//        return sb.toString();
//    }


}

class TagTaxaDistribution {
    private final Tag myTag;
    private final TaxaDistribution myTaxaDist;

    public TagTaxaDistribution(Tag myTag, TaxaDistribution myTaxaDist) {
        this.myTag = myTag;
        this.myTaxaDist = myTaxaDist;
    }

    public Tag tag() {
        return myTag;
    }

    public TaxaDistribution taxaDist() {
        return myTaxaDist;
    }
}

