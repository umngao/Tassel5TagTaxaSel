/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
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
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;

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
    private PluginParameter<Integer> refKmerLen = new PluginParameter.Builder<Integer>("refKmerLen", 300, Integer.class).guiName("Reference Kmer Length")
            .description("Length of kmers created from reference genome to store in the DB. This should be at as long or longer thant the kmerLen parameter used for storing input sequence tags.").build();
    
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
            if (tagsWithDepth.isEmpty()) {
                System.out.println("\nNo tags found with minimum depth " + minCount() + ". Halting Run.\n");
                return null;
            }
            
            System.out.println("Calling createKmerSeedsFromDBTags");
           // Create map of kmer seeds from db tags
            createKmerSeedsFromDBTags(tagsWithDepth, kmerTagMap,  window);
            System.out.println("TotalTime for createKmerSeedsFromDBTags was " + (System.nanoTime() - time) / 1e9 + " seconds");
 
            System.out.println("Size of tagsWithDepth: " + tagsWithDepth.size());
            System.out.println("Size of kmerTagMap keyset: " + kmerTagMap.keySet().size());
            
            // LCJ - this is for debug           
            int[] minMaxKmerTagCount = getMaxValueAtKey(kmerTagMap);
            System.out.println("Max number of tags at a kmer: " + minMaxKmerTagCount[1] + ", Min number of tags at kmer: " + minMaxKmerTagCount[0]);
            
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
            
            // Make array of bitmaps for the chromosomes:  OpenBitSet polybits = new OpenBitSet(nsites);
            OpenBitSet[] chromBitMaps = new OpenBitSet[10];
            System.out.println("LCJ - start making array of chrom bitmaps ");
            for (Chromosome chrom : chromsInRef){
                
                // LCJ - remove this code - is for skipping scaffolds when using zea maize full ref genome
                if (chrom.getChromosomeNumber() == Integer.MAX_VALUE) {
                    System.out.println("RepGenAligner:processData: skipping ref genome locus: " + chrom.getName());
                    continue; // skipping scaffolds
                }               
                // REMEMBER !! If only chrom 9 is set, also change to chrom 9 below.               
               // if (chrom.getChromosomeNumber() != 9) continue; // just for initial testing !!! - remove
                // LCJ - end skip code
                
                int kmersForChrom = 0;
                int chromLength = myRefSequence.chromosomeSize(chrom);
                System.out.println("\nChecking reference chrom " + chrom.getName() + ", size: " + chromLength + " for tag kmer matches.");
                
                // create bitmap for this chrom
                OpenBitSet chromBits = new OpenBitSet(chromLength);
                
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
                                               
                        // set position in the bitmap
                        chromBits.fastSet(chromIdx);
                    }
                    // after processing, slide up by 1                                  
                    chromIdx++;
                }
                chromBitMaps[chrom.getChromosomeNumber()-1] = chromBits;
                System.out.println("Total tag seeds matching to kmers in chrom " + chrom.getName() + ": " 
                    + kmersForChrom + ", total fastBits set via cardinality: " + chromBits.cardinality());
            }
            // The bitmaps of positions matching a kmer seed start have been calculated for each chrom.
            // Take these maps and find where hits are clustered.  Determine a start positions for these
            // and attempt to align the tags.  We are aligning a matrix of full tags to chrom positions.
            
            // Hashmap holding start positions, on each chromosome, where clusters
            // of hits occur, ie kmers map to these regions.  THe map is keyed by chromosome number,
            // map values the position at the center of hit clusters.  SHould this map include the tag?  
            // For now, create it
            // later.  If we have species with non-numeric chromosomes, this would need to be text as a key.
            
            System.out.println("LCJ - entering loop to find clusters" );
            Multimap<Integer,Integer> chromMaximaMap = HashMultimap.create(); // holds chrom/positions for creating refTag
            for (int mapidx = 8; mapidx < 9; mapidx++) { // use when testing only chrom 9 !!  REPLACE ORIGINAL
            //for (int mapidx = 0; mapidx < 10; mapidx++) {
                int chrom = mapidx+1; // chroms 1-10 stored in maps 0-9
                OpenBitSet chromHits = chromBitMaps[mapidx]; // chrom 1 stored in position 0
                long chromSize = chromHits.size(); //
                
                // This holds the moving sums for specific positions.  The key is
                // the sum amount.  Only sums that are equal or greater than the minCount()
                // parameter are stored.  The values are all positions that have value X
                // as their moving-sum value.
                Multimap<Short,Integer> movingSumMap = HashMultimap.create();
 
                System.out.println("LCJ - size of OpenBiSet for chrom " + chrom + ": " 
                + chromSize + ", cardiality is: " + chromHits.cardinality());
 
                // loop through the chromoseom.  Start position is ref
                int rangeSize = refKmerLen();               
                // THis is hte position where we store the first kmer
                int firstRange = rangeSize/2; // refKmerLen shoudl be even number, but is ok if it isn't                
                short firstRangeCount = calculateFirstRangeCount(rangeSize,chromHits);
                System.out.println("LCJ - firstRangeCount: " + firstRangeCount);
                
                short currentTotal = firstRangeCount;
                
                // If the range is 10, and we have 20 values.  The last range has indices 10-19,
                // and we want to store in position 15.  When getting the last value to
                // add, it is 15 + 5 -1 = 19;
                int lastRange = (int)chromSize - firstRange ;
                // calculate the rest of the values
                for (int idx = firstRange+1; idx < lastRange; idx++) {
                    // we move up one.  If the bitMap was set for the position
                    // we drop, then drop one from the count.  If the bitMap
                    // is set for the position we add, add 1 to the count.
                    // if our range is 10, the first range is 0-9
                    // We want to store the value in the middle, which is range/2 = 5
                    // That was done above.  Here, we start writing to position 6.
                    // We want the totals now in range 1-10.  We dropped position 0,
                    // which is 6 - (10/2) -1 = 6-5-1 = 0.
                    // We want to add position 10, which is 6 + (10/2) -1,
                    // which is 6 + 5 - 1 = 10
                    int posToDrop = idx - (rangeSize/2) - 1;
                    int posToAdd = idx + (rangeSize/2) -1;
                    if (chromHits.fastGet(posToDrop)) currentTotal--;
                    if (chromHits.fastGet(posToAdd)) currentTotal++;

                    // Only store positions with hit counts where the hit count meets our minimum
                    if (currentTotal >= minCount()) {
                        movingSumMap.put(currentTotal, idx);
                    }
                }   
                
                // Now find the maximas
                // Get the keys from movingSumMap.  Sort them
                // FOr each key on the list, get their values.  Sort them.
                // To find the maxima:
                //   FOr each set of sorted positions (map values) having a specific moving count:
                //       Find the start and end of each string of consecutive positions
                //       Find the mid-point of that start/end range.
                //       Store the mid-point position as a maxima
                //   Additional checking to help weed out duplicates (eg. when moving ave drops from 7 to 6)
                //   Don't count positions where there are less than 4 consecutive positions with same value
                //   Will this make a difference?
                
                Set<Short> movingSums = movingSumMap.keySet();
                
                List<Short> sumList = new ArrayList<Short>(movingSums); // does this need to be sorted?
                sumList.stream().forEach(sum-> {
                    Collection<Integer> values = movingSumMap.get(sum);
                    List<Integer> positionsWithSum = new ArrayList<Integer>(values);
                    Collections.sort(positionsWithSum);
                    System.out.println("LCJ - total number of positions with minimum count: " + sum + " is "+ positionsWithSum.size());
                    storeToChromMaximaMap( positionsWithSum, chromMaximaMap,  chrom);
                });
                
                // write out to a file for debug
                String outFile = "/Users/lcj34/notes_files/repgen/junit_out/chrom" + chrom + "hitMaximaPositions.txt";
                BufferedWriter bw = Utils.getBufferedWriter(outFile);
                // Print the positions for this chrom
                Collection<Integer> maximaForChrom = chromMaximaMap.get(chrom);
                List<Integer> chromValues = new ArrayList<Integer>(maximaForChrom);
                Collections.sort(chromValues);
                System.out.println("LCJ - total number of maxima positions for chrom " + chrom 
                        + " is " + maximaForChrom.size() + ". Writing them to out file .... ");
                
                try {
                    StringBuilder sb = new StringBuilder();
                    int count = 0;
                    for (int idx = 0; idx < chromValues.size(); idx++) {
                        String hits = Integer.toString(chromValues.get(idx));
                        sb.append(hits);
                        sb.append("\n");
                        count++;
                        if (count > 100000) {
                            bw.write(sb.toString());
                            count = 0;
                            sb.setLength(0);
                        }                       
                    }
                    if (count > 0) {
                        // write last lines
                        bw.write(sb.toString());
                    }
                    bw.close();
                } catch (IOException ioe) {
                    System.out.println("LCJ - exception writing file " + outFile);
                    ioe.printStackTrace();
                } 
            }
            
            // This is just for debug
           // for (int mapidx = 0; mapidx < 10; mapidx++){
            for (int mapidx = 8; mapidx < 9; mapidx++) {
                int numPosForChrom = chromMaximaMap.get(mapidx+1).size();
                int chrom = mapidx+1;
                System.out.println("LCJ - number of clustered hits for chrom " + chrom + " is " + numPosForChrom);
            }
            
            // TODO:  Now that we have the positions, get the tags and then run SW 
            // Nothing below here is finished .
            
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
            
            ((RepGenSQLite)repGenData).close();  //todo autocloseable should do this but it is not working.

            myLogger.info("Finished RepGenAlignerPlugin\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }
    
    private int[] getMaxValueAtKey(Multimap<String,Tag>kmerTagMap){
       int[] minMaxTagCounts = {10000,0};
       for(String key: kmerTagMap.keySet()) {
           int count = kmerTagMap.get(key).size();
           if (count > minMaxTagCounts[1]) minMaxTagCounts[1] = count;
           if (count < minMaxTagCounts[0]) minMaxTagCounts[0] = count;
       }
       return minMaxTagCounts;
    }
    
    private short calculateFirstRangeCount(int rangeSize, OpenBitSet obs) {
        short totalhits = 0;
        for (int idx = 0; idx < rangeSize; idx++) {
            if (obs.get(idx)) totalhits++;
        }
        return totalhits;
    }
    
    private void storeToChromMaximaMap(List<Integer> positionsWithSum,Multimap<Integer,Integer> chromMaximaMap, int chrom) {
        int start = positionsWithSum.get(0);
        int previous = start;
        int next = start;
        // Is using lastRange as above correct here?
        for (int idx = 1; idx < positionsWithSum.size();) {
            next = positionsWithSum.get(idx);
            while (next == previous+1){
                previous = next;
                idx++;
                if (idx < positionsWithSum.size()) next = positionsWithSum.get(idx);
                else break;
            }
            // Hoping to skip some of the overlap positions, e.g. when the sum at position X is 5
            // and the sum at position X+1 is 4, I don't want both X and X+1 in the file
            // Having this code dropped the number of entries from 12830 to 4258
            int haveEnough = (previous-start)/2;
            if (haveEnough >= 2) { // requires 4 positions with same hit count in a row.
                // Get the the middle of the positions, then add 8.  This is because the
                // kmer is 16bp long and we want to center in the middle of it
                int maxima = start + haveEnough + 8;  
                chromMaximaMap.put(chrom, maxima);
            }
            start = next;
            previous = next;
            idx++;
        }
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
    /**
     * Length of  kmers as tag sequences in the db
     *
     * @return kmerLen
     */
    public Integer refKmerLen() {
        return refKmerLen.value();
    }

    /**
     * Set Kmer length. Length of kmers to be stored
     * as tag sequences in the db
     *
     * @param value kmer length
     *
     * @return this plugin
     */
    public RepGenAlignerPlugin refKmerLen(Integer value) {
        refKmerLen = new PluginParameter<>(refKmerLen, value);
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

