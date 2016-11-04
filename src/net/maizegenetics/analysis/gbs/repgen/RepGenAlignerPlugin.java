/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;
//import org.biojava.nbio.alignment.SmithWaterman;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.SmithWaterman;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;

/**
 * This plugin takes an existing repGen db, grabs the tags
 * whose depth meets that specified in the minCount parameter,
 * and makes a best guess at each tag's alignment against a
 * reference genome.
 * 
 * kmer seeds are created from the stored tags.  The ref genome
 * is walked with a sliding window of 1.  Smith Waterman is used
 * to determine alignment score.
 * 
 * @author lcj34
 *
 */
public class RepGenAlignerPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(RepGenAlignerPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<String> refGenome = new PluginParameter.Builder<String>("ref", null, String.class).guiName("Reference Genome File").required(true)
            .description("Referemce Genome File for aligning against ").build();
    private PluginParameter<Integer> myMinCount = new PluginParameter.Builder<Integer>("minCount", 1, Integer.class).guiName("Min Count")
            .description("Minimum count of reads for a tag to be aligned").build();
    private PluginParameter<Integer> seedLen = new PluginParameter.Builder<Integer>("seedLen", 16, Integer.class).guiName("Seed Kmer Length")
            .description("Length of kmer seed created from DB tags and used as seeds for aligning against the reference genome.").build();
    private PluginParameter<Integer> kmerLen = new PluginParameter.Builder<Integer>("kmerLen", 150, Integer.class).guiName("Kmer Length")
            .description("Length of kmer from fastq reads stored as the tag sequence in the DB.").build();
    
    static GenomeSequence myRefSequence = null;

    public RepGenAlignerPlugin() {
        super(null, false);
    }

    public RepGenAlignerPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RepGenAlignerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenAlignerPlugin: postProcessParameters: Input DB not set or found");
        }
        if (!refGenome.isEmpty()) {
            myRefSequence = GenomeSequenceBuilder.instance(refGenome());
        } else {
            throw new IllegalArgumentException("RepGenAlignerPlugin: postProcessParameters: reference genome not set or found");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
        // Make seeds from the tags with minimum counts in the db
        // Keep map of kmer/tags, align the kmers to the reference genome
        // if kmer seed matches, attempt alignment of tag at this position
        // store best alignment for each tag, load data into TagMapping and
        // PhysicalMapPosition
        try {           
            System.out.println("RepGenAlignerPlugin:processData begin"); 
            RepGenDataWriter repGenData=new RepGenSQLite(inputDB());
                       
            Multimap<String,Tag> kmerTagMap = HashMultimap.create();
            int window = 20;
            
            Map<Tag, Integer> tagsWithDepth = repGenData.getTagsWithDepth(minCount());
            
            System.out.println("Calling createKmerSeedsFromDBTags");
           // Create map of kmer seeds from db tags
            createKmerSeedsFromDBTags(tagsWithDepth, kmerTagMap,  window);
            System.out.printf("TotalTime for createKmerSeedsFromDBTags was %g sec%n", (double) (System.nanoTime() - time) / 1e9);
            
            time = System.nanoTime();
            // get reference genome
            myRefSequence = GenomeSequenceBuilder.instance(refGenome());
            
            Multimap<Tag,AlignmentInfo> tagAlignInfoMap = HashMultimap.create();
            // For each chrom in refSequence, walk the refSequence looking for kmers
            // matching those in the kmerTagMap.  When kmer found, align the tags for that
            // kmer against the reference.  Of each tag's alignment choices, store the best
            // one in the TagMapping table
            Set<Chromosome> chromsInRef = myRefSequence.chromosomes();
            int kmersFound = 0;
            for (Chromosome chrom : chromsInRef){
                // LCJ - remove this code - is for skipping scaffolds when using zea maize full ref genome
                if (chrom.getChromosomeNumber() == Integer.MAX_VALUE) continue; // skipping scaffolds
                // LCJ - end skip code
                int kmersForChrom = 0;
                int chromLength = myRefSequence.chromosomeSize(chrom);
                System.out.println("\nChecking reference chrom " + chrom.getName() + " for tag kmer matches.");
                time = System.nanoTime();
                for (int chromIdx = 0; chromIdx < chromLength;) {
                                        
                    // chromosomeSequence:  start and end are inclusive and 1-based, adjust for this in call below
                    int end = Math.min(chromLength, chromIdx+seedLen());
                    byte[] chromSeedBytes = myRefSequence.chromosomeSequence(chrom,chromIdx+1, end);
                    // we're not processing from fastq, we're processing from tag table in db
                    int badValue = checkForN(chromSeedBytes);
                    if (badValue >=0) {
                        // found a non-ACGT character,  re-set index and
                        // grab new kmer from the chromosome
                        chromIdx += badValue+1; // returned value was at bad base - move beyond it
                        continue;
                    }
                                       
                    // Convert back to allele string
                    String chromKmerString = NucleotideAlignmentConstants.nucleotideBytetoString(chromSeedBytes);
                    // chromKmerBytes were good - check if this kmer exists among our seeds
                    if (kmerTagMap.containsKey(chromKmerString)){
                        kmersFound++;
                        kmersForChrom++;
                        // get full size of tag - this returns it in numeric byte form
                        byte[] chromTagBytes = myRefSequence.chromosomeSequence(chrom,chromIdx+1,chromIdx + kmerLen());
                        String chromTagAlleles = NucleotideAlignmentConstants.nucleotideBytetoString(chromTagBytes);
                        
                        // for all tags containing this kmer, run SMith Waterman to align.
                        // return values in TagAlignmentMap                       
                        Collection<Tag> tagsToAlign = kmerTagMap.get(chromKmerString);
                        List tagsList = new ArrayList(tagsToAlign);
                        calculateAlignmentInfo(chromTagAlleles.getBytes(),tagsList,tagAlignInfoMap,chrom.getName(),chromIdx+1); // position is 1 based
                    }
                    // after processing, slide up by 1                  
                    time = System.nanoTime();
                    chromIdx++;
                }
                System.out.println("Aligning to reference chrom " + chrom.getName() + " took " + (System.nanoTime() - time)/1e9 + " seconds");
                System.out.println("Total tag seeds matching to kemrs in chrom " + chrom.getName() + ": " + kmersForChrom);
            }
            
            // map of tag/positions to store in the db
            // Ed's notes show a multimap, but I thought each tag got 1 position best
            // on the phsical map.  Here, we run through  each tag's ALignmentInfo data and
            // store the one with the best score.
            Multimap<Tag,Position> kmerPositionMap =  HashMultimap.create();
           
            Set<Tag> tagSet= tagAlignInfoMap.keySet();
            Iterator<Tag> tagIterator = tagSet.iterator();
            
            System.out.println("\nAlignments are finished, total number of kmers matched in all chromosomes:" + kmersFound);
            System.out.println("\n Start loop to select best alignments");
            while (tagIterator.hasNext()) {
                Tag tag = (Tag)tagIterator.next();
                
                List<AlignmentInfo>  alignments = (List<AlignmentInfo>) tagAlignInfoMap.get(tag);
                
                AlignmentInfo bestAlignment = alignments.get(0);
                int bestScore = -1;
                for (AlignmentInfo ai : alignments) {
                    if (ai.score() > bestScore) {
                        bestAlignment = ai;
                        bestScore = bestAlignment.score();
                    }
                }
                Chromosome chromosome = new Chromosome(bestAlignment.chromosome());
                Position position=new GeneralPosition
                        .Builder(chromosome,bestAlignment.position())
                        .strand((byte)1) // default to forward for now
                        //.strand((byte)1)
                        .addAnno("mappingapproach", "SmithWaterman")
                        .build();
                // create Position object from alignment data
                kmerPositionMap.put(tag, position);
            }
            
            System.out.println("\n Adding data to repGenDB");
            // add reference and mapping approach to db
            repGenData.addReferenceGenome(refGenome());
            repGenData.addMappingApproach("SmithWaterman");
            // Add tags and alignment to db
            repGenData.putTagAlignments(kmerPositionMap, refGenome());
            
            ((TagDataSQLite)repGenData).close();  //todo autocloseable should do this but it is not working.

            myLogger.info("Finished RepGenAlignerPlugin\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }
    
    private void createKmerSeedsFromDBTags(Map<Tag, Integer> tagWithDepth,
            Multimap<String,Tag> kmerTagMap, int window) {
        for (Tag tag : tagWithDepth.keySet()) {
            String tagSequence = tag.sequence();
            int maxIdx = tagSequence.length() - window;
            for (int seqIdx = 0; seqIdx < maxIdx;) {
                String kmer = tagSequence.substring(seqIdx, seqIdx + seedLen());  
                byte[] kmerBytesAsNum = NucleotideAlignmentConstants.convertHaplotypeStringToAlleleByteArray(kmer);
                
                int badValues = checkForN(kmerBytesAsNum);
                if (badValues >=0) {
                    // found a non-ACGT character,  re-set index and
                    // grab new kmer
                    seqIdx += badValues+1;
                    continue;
                }
                // good values - add unique kmer to map
                kmerTagMap.put(kmer, tag);
                // add reverse complement of kmer
                // THis isn't right ... this tag does not have this KMER.
                // And we need to know it is a reverse complement.  Ed's code
                // creates rc from the ref genome, and stores the position as negative
                // We are creating seeds from tags, not the ref genomes.  how to include RC of kmer?
                byte[] kmerRC = NucleotideAlignmentConstants.reverseComplementAlleleByteArray(kmer.getBytes());
                String kmerRCString = new String(kmerRC);
                kmerTagMap.put(kmerRCString, tag); 
                // loop end - increment seqIdx for next kmer
                seqIdx += window;
            }               
        }
    }
    
    private void calculateAlignmentInfo(byte[] refTagSequence,List<Tag> tags, 
            Multimap<Tag,AlignmentInfo> tagAlignInfoMap, String chrom, int position){
        // For each tag on the tags list, run SW against it and the refTagSequence
        // Store alignment info in the map.
        
        for (Tag tagToMatch : tags) {
            SmithWaterman mySM = new SmithWaterman(refTagSequence, tagToMatch.sequence().getBytes());
            int score = mySM.computeScore();
            AlignmentInfo tagAI = new AlignmentInfo(Arrays.toString(refTagSequence),chrom,position,score);
            tagAlignInfoMap.put(tagToMatch,tagAI);
        }       
    }
    
    // Using edit distance as an alignment score may be faster than SW,
    // which is taking way too long (ie, not finishing) when run from my laptop.
    private int computeEditDistance(String ncTagSeq, String cTagSeq) {
        int cTagLen = cTagSeq.length();
        int ncTagLen = ncTagSeq.length();
        boolean cTagIsLongest = cTagLen > ncTagLen? true : false;
        int maxCompare = cTagIsLongest ? ncTagLen : cTagLen;

        int editDistance1 = 0;
        int editPosition = 0;
        for ( int idx1 = 0; idx1 < maxCompare; idx1++) {
            if (ncTagSeq.charAt(idx1) != cTagSeq.charAt(idx1))
                editDistance1++;
            if (editDistance1== 1) editPosition = idx1; // save 1st mismatch position
        }
        if (editDistance1 == 1) return editDistance1; // 1 mismatch is best possible

        // At first mismatched bp, shift the non-canonical tag left 1 position and re-align
        // This checks for a 1 bp insertion
        int editDistance2 = 1;  // accounts for the gap
        for (int idx2 = editPosition; idx2 < maxCompare-1; idx2++ ) { // account for difference at beginning (handle 0 to negative)
            if (ncTagSeq.charAt(idx2+1) != cTagSeq.charAt(idx2))
                editDistance2++;
        }      
        if (editDistance2 == 1) return editDistance2;
        int bestEditDistance = editDistance1 < editDistance2 ? editDistance1 : editDistance2;

        // At first mismatch bp, Shift the non-canonical tag right 1 position and re-align
        // This checks for a 1 bp deletion
        int editDistance3 = 1;  // accounts for the gap
        for (int idx3 = editPosition; idx3 < maxCompare-1; idx3++ ) { //
            if (ncTagSeq.charAt(idx3) != cTagSeq.charAt(idx3+1))
                editDistance3++;
        }

        // Compute which is shortest, editDistance1, 2 or 3 and return that value
        return bestEditDistance < editDistance3 ? bestEditDistance : editDistance3;
    }
    // Checks for non ACGT character in kmer.  If found, returns the deviate
    // position in the array.  If not found, -1 is returned. -1 is good return
    private int checkForN(byte[] kmer) {
        boolean flag = false;
        int idx;
       // byte[] kmerInBytes = NucleotideAlignmentConstants.convertHaplotypeStringToAlleleByteArray(kmer);
        for (idx = 0; idx < kmer.length; idx++) {
            if (kmer[idx] >3) { // non ACGT value
                flag = true;
                break;
            }
        }
        if (flag) return idx; // toss sequences with non ACGT value
        else return -1;
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
    public RepGenAlignerPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Output fastq file to use as input for BWA or bowtie2
     *
     * @return Output File
     */
    public String refGenome() {
        return refGenome.value();
    }

    /**
     * Set Output File. Output fastq file to use as input
     * for BWA or bowtie2
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin refGenome(String value) {
        refGenome = new PluginParameter<>(refGenome, value);
        return this;
    }

    /**
     * Minimum count of reads for a tag to be output
     *
     * @return Min Count
     */
    public Integer minCount() {
        return myMinCount.value();
    }

    /**
     * Set Min Count. Minimum count of reads for a tag to
     * be output
     *
     * @param value Min Count
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin minCount(Integer value) {
        myMinCount = new PluginParameter<>(myMinCount, value);
        return this;
    }
    
    /**
     * Length of seed kmers
     *
     * @return seed len
     */
    public Integer seedLen() {
        return seedLen.value();
    }

    /**
     * Set Seed kmer length. Length of seeds to 
     * use in aligning.
     *
     * @param value seed length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin seedLen(Integer value) {
        seedLen = new PluginParameter<>(seedLen, value);
        return this;
    }
    /**
     * Length of  kmers as tag sequences in the db
     *
     * @return kmerLen
     */
    public Integer kmerLen() {
        return kmerLen.value();
    }

    /**
     * Set Kmer length. Length of kmers to be stored
     * as tag sequences in the db
     *
     * @param value kmer length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin kmerLen(Integer value) {
        kmerLen = new PluginParameter<>(kmerLen, value);
        return this;
    }
}

class AlignmentInfo implements Comparable<AlignmentInfo>{
    private final String refSequence;
    private  final String myChromosome;
    private  final int myPosition;
    private  final int myScore;

    public AlignmentInfo(String refSequence, String chromosome, int position, int score) {
        this.refSequence = refSequence;
        this.myChromosome = chromosome;
        this.myPosition = position;
        this.myScore = score;
    }

    public String chromosome() {
        return myChromosome;
    }

    public int position() {
        return myPosition;
    }
    
    public  int score() {
        return myScore;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Alignment:");
        sb.append("\tRefSequence:").append(refSequence);
        sb.append("\tChr:").append(myChromosome);
        sb.append("\tPos:").append(myPosition);
        sb.append("\tName:").append(myScore);
        sb.append("\n");
        return sb.toString();
    }
    @Override
    public int compareTo(AlignmentInfo other) {
        // TODO Auto-generated method stub
        return this.myScore > other.score() ? 1: this.myScore < other.score() ? -1 : 0;
    }
}

