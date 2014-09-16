/*
 * DiscoverySNPCallerPluginV2
 */
package net.maizegenetics.analysis.gbs.v2;

import com.google.common.collect.*;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.tag.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.util.ConcurrencyTools;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import net.maizegenetics.plugindef.PluginParameter;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.*;

/**
 * This class aligns tags at the same physical location against one another,
 * calls SNPs, and then outputs the SNPs to a HapMap file.
 *
 * It is multi-threaded, as there are substantial speed increases with it.
 *
 * @author Ed Buckler
 * @author Jeff Glaubitz
 */
public class DiscoverySNPCallerPluginV2 extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(DiscoverySNPCallerPluginV2.class);

    private PluginParameter<String> myInputDB = new PluginParameter.Builder<>("i", null, String.class).guiName("Input GBS Database").required(true).inFile()
            .description("Input Database file if using SQLite").build();
    private PluginParameter<Double> myMinMinorAlleleFreq = new PluginParameter.Builder<>("mnMAF", 0.01, Double.class).guiName("Min Minor Allele Freq")
            .description("Minimum minor allele frequency").build();
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

    private TagDataWriter tagDataWriter = null;
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
        myLogger.info("Finding SNPs in " + inputDB() + ".");
        myLogger.info(String.format("StartChr:%d EndChr:%d %n", startChromosome(), endChromosome()));
        //DataOutputStream locusLogDOS = openLocusLog(logFile());
        if (customSNPLogging) {
//            myCustomSNPLog = new CustomSNPLog(logFile());
        }
        for (int chr = startChromosome(); chr <= endChromosome(); chr++) {
            myLogger.info("\n\nProcessing chromosome " + chr + "...");
            if (includeReference) {
                //refGenomeChr = readReferenceGenomeChr(referenceGenomeFile(), chr);
                if (refGenomeChr == null) {
                    myLogger.info("  WARNING: chromosome " + chr + " not found in the reference genome file. Skipping this chromosome.");
                    continue;
                }
            }
            discoverSNPsOnChromosome(chr);
            myLogger.info("Finished processing chromosome " + chr + "\n\n");
        }
        ConcurrencyTools.shutdown();
        return null;
    }

    @Override
    public void postProcessParameters() {

        if (myInputDB.isEmpty()) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Input Tags by Taxa File not Set.");
        } else {
            tagDataWriter =new TagDataSQLite(inputDB());
        }
        if (tagDataWriter == null) {
            throw new IllegalArgumentException("DiscoverySNPCallerPlugin: postProcessParameters: Problem reading Tags by Taxa File: " + inputDB());
        }
        if (!myRefGenome.isEmpty()) {
            includeReference = true;
        }

        if (callBiallelicSNPsWithGap() && includeGaps()) {
            throw new IllegalArgumentException("The callBiSNPsWGap option is mutually exclusive with the inclGaps option.");
        }

        if (endChromosome() - startChromosome() < 0) {
            throw new IllegalArgumentException("The start chromosome is larger than the end chromosome.");
        }

        myLogger.info(String.format("MinMAF:%g %n", minMinorAlleleFreq()));
        myLogger.info(String.format("includeRare:%s includeGaps:%s %n", includeRareAlleles(), includeGaps()));
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

    public void discoverSNPsOnChromosome(int targetChromo) {
        int siteCnt = 0;
        tagDataWriter.getCutPositionTagTaxaMap(new Chromosome("" + targetChromo), -1, -1).entrySet().stream()
                .forEach(emp -> findAlleleByAlignment(emp.getKey(),emp.getValue()));
        myLogger.info("Number of marker sites recorded for chr" + targetChromo + ": " + siteCnt);
    }

    Table<Position, Byte, List<TagTaxaDistribution>> findAlleleByAlignment(Position cutPosition, Map<Tag,TaxaDistribution> tagTaxaMap) {
        System.out.println("cutPosition = [" + cutPosition + "], tagTaxaMap = [" + tagTaxaMap + "]");
        if(tagTaxaMap.isEmpty()) return null;  //todo why would this be empty?
        final int numberOfTaxa=tagTaxaMap.values().stream().findFirst().get().maxTaxa();
        if(tagTaxaMap.size()<2) {//homozygous
            return null;  //consider reporting homozygous
        }
        if(!tagTaxaMap.keySet().stream().anyMatch(Tag::isReference)) {
            tagTaxaMap=setCommonToReference(tagTaxaMap);
        }
        final double taxaCoverage=tagTaxaMap.values().stream().mapToInt(TaxaDistribution::numberOfTaxaWithTag).sum()/(double)numberOfTaxa;  //todo this could be changed to taxa with tag
        if(taxaCoverage < myMinLocusCoverage.value()) {
            return null;  //consider reporting low coverage
        }
        Map<Tag,String> alignedTags=alignTags(tagTaxaMap.keySet());
        Table<Position, Byte, List<TagTaxaDistribution>> tAlign=convertAlignmentToTagTable(alignedTags, tagTaxaMap,  cutPosition);
        List<Position> positionToKeep=tAlign.rowMap().entrySet().stream()
                .filter(entry -> (double) numberTaxaAtSiteIgnoreGaps(entry.getValue()) / (double) numberOfTaxa > minLocusCoverage())
                .filter(entry -> {
                    List<Tuple<Byte,Integer>> aC=alleleTaxaCounts(entry.getValue());
                    return (aC.size()>1 && aC.get(1).y>(minMinorAlleleFreq()*taxaCoverage));
                })
                .map(Map.Entry::getKey) //get Position
                .collect(Collectors.toList());
        //todo convert to stream
        Multimap<Tag,Allele> tagAllelemap= HashMultimap.create();
        for (Position position: positionToKeep) {
            for (Map.Entry<Byte, List<TagTaxaDistribution>> entry : tAlign.row(position).entrySet()) {
                Allele allele=new SimpleAllele(entry.getKey(),position);
                for (TagTaxaDistribution tagTaxaDistribution : entry.getValue()) {
                    tagAllelemap.put(tagTaxaDistribution.tag(),allele);
                }
            }
        }
        tagDataWriter.putTagAlleles(tagAllelemap);
        return tAlign;
    }

    private static Map<Tag,TaxaDistribution> setCommonToReference(Map<Tag,TaxaDistribution> tagTaxaMap) {
        Tag commonTag=tagTaxaMap.entrySet().stream()
                .max(Comparator.comparingInt(e -> e.getValue().numberOfTaxaWithTag()))
                .map(e -> e.getKey())
                .get();
        TaxaDistribution commonTD=tagTaxaMap.get(commonTag);
        Tag refTag=TagBuilder.instance(commonTag.seq2Bit(),commonTag.seqLength()).reference().build();
        ImmutableMap.Builder<Tag, TaxaDistribution> tagTaxaMapBuilder=new ImmutableMap.Builder<>();
        for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaMap.entrySet()) {
            if(entry!=commonTag) {tagTaxaMapBuilder.put(entry);}
            else {tagTaxaMapBuilder.put(refTag,commonTD);}
        }
        return tagTaxaMapBuilder.build();
    }

    private static List<Tuple<Byte,Integer>> alleleTaxaCounts(Map<Byte, List<TagTaxaDistribution>> alleleDistMap) {
        List<Tuple<Byte,Integer>> alleleCnt = new ArrayList<>();
        for (Map.Entry<Byte, List<TagTaxaDistribution>> entry : alleleDistMap.entrySet()) {
            alleleCnt.add(new Tuple<>(entry.getKey(),numberOfTaxa(entry.getValue())));
        }
        Collections.sort(alleleCnt, (Tuple<Byte, Integer> a, Tuple<Byte, Integer> b) -> b.y.compareTo(a.y));  //a,b reversed so common allele  first
        return alleleCnt;
    }

    /**
     * Aligns a set of tags anchored to the same reference position.
     * @return map with tag(values) mapping to String with alignment
     */
    static Map<Tag,String> alignTags(Collection<Tag> tags) {
        List<DNASequence> lst;
        lst = tags.stream()
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
     * Converts at TagAlignment to a Guava Table with
     * @param alignedTags
     * @param tagTaxaDistMap
     * @param refStartPosition
     * @return
     */
    static Table<Position, Byte, List<TagTaxaDistribution>> convertAlignmentToTagTable(Map<Tag,String> alignedTags,
                        Map<Tag,TaxaDistribution> tagTaxaDistMap, Position refStartPosition) {
        Table<Position, Byte, List<TagTaxaDistribution>> alignT= TreeBasedTable.create(); //These could be sorted by depth
        final List<Optional<Position>> referencePositions=referencePositions(refStartPosition,alignedTags);
        alignedTags.forEach((t,s) -> {
            TagTaxaDistribution td=new TagTaxaDistribution(t,tagTaxaDistMap.get(t));
            System.out.println("convertAlignmentToTagTable:"+td);
            for (int i = 0; i < s.length(); i++) {
                byte allele=getNucleotideAlleleByte(s.charAt(i));
                Optional<Position> p=referencePositions.get(i);
                if(!p.isPresent()) continue;  //skip alignment sites not in reference
                List tdList=alignT.get(p.get(),allele);
                if(tdList!=null) {tdList.add(td);}
                else {
                    List<TagTaxaDistribution> ttdL=new ArrayList<>();
                    ttdL.add(td);
                    alignT.put(p.get(),allele, ttdL);
                }
            }
        });
        return alignT;
    }

    public static String toString(Table<Position, Byte, List<TagTaxaDistribution>> snpTagTable) {
        StringBuilder sb=new StringBuilder();
        sb.append(String.format("Rows %d Columns %d\n", snpTagTable.rowKeySet().size(), snpTagTable.columnKeySet().size()));
        sb.append("Columns bases:");
        snpTagTable.columnKeySet().stream().forEach(b -> sb.append(b+","));
        sb.append("\n");
        snpTagTable.rowMap().forEach((p, ttd) -> {
            sb.append(String.format("Position: %d", p.getPosition()));
            for (Map.Entry<Byte, List<TagTaxaDistribution>> byteListEntry : ttd.entrySet()) {
                sb.append(String.format("\tAllele: %d ->", byteListEntry.getKey()));
                byteListEntry.getValue().forEach(tt -> sb.append(String.join(",",tt.tag().sequence()+"#"+tt.taxaDist().numberOfTaxaWithTag()+" ")));//not doing this correctly
            }
            sb.append("\n");

        });
        return sb.toString();
    }

    private static int numberOfTaxa(List<TagTaxaDistribution> tagTaxaDistributions) {
        return tagTaxaDistributions.stream().mapToInt(t -> t.taxaDist().numberOfTaxaWithTag()).sum();
    }

    private static int numberTaxaAtSite(Map<Byte,List<TagTaxaDistribution>> tagTaxaDistributions) {
        return tagTaxaDistributions.values().stream()
                .flatMap(tList -> tList.stream())
                .mapToInt(t -> t.taxaDist().numberOfTaxaWithTag())
                .sum();
    }

    private static int numberTaxaAtSiteIgnoreGaps(Map<Byte,List<TagTaxaDistribution>> tagTaxaDistributions) {
        return tagTaxaDistributions.entrySet().stream()
                .filter(e -> e.getKey()!=GAP_ALLELE)
                .flatMap(tList -> tList.getValue().stream())
                .mapToInt(tag -> tag.taxaDist().numberOfTaxaWithTag())
                .sum();
    }

//    private void removePosition(Table<Position, Byte, List<TagTaxaDistribution>> table, Position posToRemove){
//        table.row(posToRemove).forEach((b,ttd) -> table.remove(posToRemove,b));
//    }

    private static List<Optional<Position>> referencePositions(Position refStartPosition, Map<Tag,String> alignedTags){
        Tag refTag=alignedTags.keySet().stream()
                .filter(Tag::isReference)
                .findFirst().orElseThrow(() -> new IllegalStateException("Reference not found"));
        AtomicInteger start=new AtomicInteger(refStartPosition.getPosition());
        return alignedTags.get(refTag).chars()
                .mapToObj(c -> (c == GAP_ALLELE_CHAR) ? Optional.<Position>empty() :
                        Optional.<Position>of(new GeneralPosition.Builder(refStartPosition.getChromosome(), start.getAndIncrement()).build()))
                .collect(Collectors.toList());
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(DiscoverySNPCallerPluginV2.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
//    // TODO: Replace <Type> with specific type.
//    public <TagDataWriter> runPlugin(DataSet input) {
//        return (TagDataWriter) performFunction(input).getData(0).getData();
//    }

    /**
     * Input TagsByTaxa file (if hdf5 format, use .hdf or
     * .h5 extension)
     *
     * @return Input Tags by Taxa File
     */
    public String inputDB() {
        return myInputDB.value();
    }

    /**
     * Set Input Tags by Taxa File. Input TagsByTaxa file
     * (if hdf5 format, use .hdf or .h5 extension)
     *
     * @param value Input Tags by Taxa File
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 inputDB(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
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
     * Minimum locus coverage (proportion of Taxa with a genotype)
     *
     * @return Min Locus Coverage
     */
    public Double minLocusCoverage() {
        return myMinLocusCoverage.value();
    }

    /**
     * Set Min Locus Coverage. Minimum locus coverage (proportion
     * of Taxa with a genotype)
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
     * Path to reference genome in fasta format. Ensures that
     * a tag from the reference genome is always included
     * when the tags at a locus are aligned against each other
     * to call SNPs. The reference allele for each site is
     * then provided in the output HapMap files, under the
     * taxon name "REFERENCE_GENOME" (first taxon). DEFAULT:
     * Don't use reference genome.
     *
     * @return Reference Genome File
     */
    public String referenceGenomeFile() {
        return myRefGenome.value();
    }

    /**
     * Set Reference Genome File. Path to reference genome
     * in fasta format. Ensures that a tag from the reference
     * genome is always included when the tags at a locus
     * are aligned against each other to call SNPs. The reference
     * allele for each site is then provided in the output
     * HapMap files, under the taxon name "REFERENCE_GENOME"
     * (first taxon). DEFAULT: Don't use reference genome.
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

    /**
     * Include the rare alleles at site (3 or 4th states)
     *
     * @return Include Rare Alleles
     */
    public Boolean includeRareAlleles() {
        return myIncludeRareAlleles.value();
    }

    /**
     * Set Include Rare Alleles. Include the rare alleles
     * at site (3 or 4th states)
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
     * Set Include Gaps. Include sites where major or minor
     * allele is a GAP
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
     * Include sites where the third allele is a GAP (mutually
     * exclusive with inclGaps)
     *
     * @return Call Biallelic SNPs with Gap
     */
    public Boolean callBiallelicSNPsWithGap() {
        return myCallBiSNPsWGap.value();
    }

    /**
     * Set Call Biallelic SNPs with Gap. Include sites where
     * the third allele is a GAP (mutually exclusive with
     * inclGaps)
     *
     * @param value Call Biallelic SNPs with Gap
     *
     * @return this plugin
     */
    public DiscoverySNPCallerPluginV2 callBiallelicSNPsWithGap(Boolean value) {
        myCallBiSNPsWGap = new PluginParameter<>(myCallBiSNPsWGap, value);
        return this;
    }



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

