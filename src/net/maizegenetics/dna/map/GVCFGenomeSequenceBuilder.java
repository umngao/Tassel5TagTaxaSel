package net.maizegenetics.dna.map;

import com.google.common.base.Splitter;
import com.google.common.collect.*;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;


/**
 * Created by zrm22 on 3/27/17.
 */
public class GVCFGenomeSequenceBuilder extends GenomeSequenceBuilder {
    private static final Pattern TAB_PATTERN = Pattern.compile("[\\t]+");

    /**
     * Builds GenomeSequence from a fasta file and a GVCF file.
     *
     * @param fastaFileName full path to fasta file
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(String fastaFileName, String gvcfFileName) throws Exception{
        Function<Character, Character> charConversion = (c) -> c;
        return instance(fastaFileName, charConversion, gvcfFileName);
    }

    /**
     * Builds GenomeSequence from a fasta file.  The char conversion provide a mechanism to convert upper and lower case
     * or convert one case to N.  This is useful if a case if used to define a certain class of bases
     *
     * @param fastaFileName  full path to fasta file
     * @param charConversion lambda Function to convert characters
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(String fastaFileName, Function<Character, Character> charConversion, String gvcfFileName) throws Exception{
        Map<Chromosome, byte[]> chromPositionMap = readReferenceGenomeChr(fastaFileName, charConversion);
        //need to create a Map<chr,Map<range,Map<AnnotationName,Value>>>
        Map<Chromosome, RangeMap<Integer,GeneralAnnotationStorage>> gvcfAnnotationsAndCalls = readGVCFFile(gvcfFileName);
        PositionList gvcfPositionsAndAnnotations = readGVCFFilePositionList(gvcfFileName);
        return new HalfByteGenomeSequenceGVCF(chromPositionMap,gvcfPositionsAndAnnotations);
    }

    private static Map<Chromosome,RangeMap<Integer,GeneralAnnotationStorage>> readGVCFFile(String gvcfFileName) {
        HashMap<Chromosome, RangeMap<Integer,GeneralAnnotationStorage>> chromosomeRangeMapHashMap = new HashMap<>();
        //TODO multithread using similar code to BuilderFromVCF

        return chromosomeRangeMapHashMap;
    }


    private static PositionList readGVCFFilePositionList(String gvcfFileName) throws Exception{
        ArrayList<Position> positionArrayList = new ArrayList<>();

        BufferedReader gvcfFileReader = new BufferedReader(new FileReader(gvcfFileName));
        //Loop through the headers
        String currentLine = "";
        while(!(currentLine = gvcfFileReader.readLine()).startsWith("#CHROM")) {

        }
        //parse the header
        String[] header = TAB_PATTERN.split(currentLine);
        HeaderPositions hp= new HeaderPositions(header);
        GVCFPositionRecord gvcfPositionRecord = new GVCFPositionRecord(hp);

        while((currentLine = gvcfFileReader.readLine())!=null) {
            positionArrayList.add(gvcfPositionRecord.parseGVCFRecords(currentLine));
        }
        PositionList instance = new PositionArrayList(positionArrayList, "AGPv3");
        return instance;
    }
}
class HeaderPositions {
    final int NUM_HAPMAP_NON_TAXA_HEADERS;
    final int GENOIDX;
    final int SNPID_INDEX;
    //  final int VARIANT_INDEX;
    final int FILTER_INDEX;
    final int QUAL_INDEX;
    final int CHROMOSOME_INDEX;
    final int POSITION_INDEX;
    final int REF_INDEX;
    final int ALT_INDEX;
    final int INFO_INDEX;
    final int FORMAT_INDEX;

    public HeaderPositions(String[] header){
        int chrIdx=firstEqualIndex(header,"#CHROM");
        if(chrIdx<0) chrIdx=firstEqualIndex(header,"#CHR");
        CHROMOSOME_INDEX=chrIdx;
        POSITION_INDEX=firstEqualIndex(header,"POS");
        SNPID_INDEX=firstEqualIndex(header,"ID");
        REF_INDEX=firstEqualIndex(header,"REF");
        ALT_INDEX=firstEqualIndex(header,"ALT");
        QUAL_INDEX=firstEqualIndex(header,"QUAL");
        FILTER_INDEX=firstEqualIndex(header,"FILTER");
        INFO_INDEX=firstEqualIndex(header,"INFO");
        FORMAT_INDEX=firstEqualIndex(header,"FORMAT");

        NUM_HAPMAP_NON_TAXA_HEADERS=Math.max(INFO_INDEX,FORMAT_INDEX)+1;
        GENOIDX=NUM_HAPMAP_NON_TAXA_HEADERS;
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i=0; i<sa.length; i++) {
            if(sa[i].equals(match)) return i;
        }
        return -1;
    }

}

/**
 * ReferenceGenomeSequence class.  This class is used to read chromosome sequences
 * from fasta files.  Data is stored as half-bytes packed into a byte array.
 * This byte array comprises the "value" for a hash map whose key is a
 * Chromosome object.
 *
 * The class also contains methods to obtain a full or partial genome sequence for a
 * specified stored chromosome.
 *
 * @author Lynn Johnson, Zack Miller
 *
 */
class HalfByteGenomeSequenceGVCF implements GVCFGenomeSequence{
    private Map<Chromosome, byte[]> chromPositionMap;
    private Map<Chromosome, Integer> chromLengthLookup=new HashMap<>();
    private RangeMap<Long,Chromosome> wholeGenomeIndexMap= TreeRangeMap.create();
    private  PositionList gvcfAnnotationsAndCalls;
    private final long genomeSize;
    private BitSet maskBitSet;
    private BitSet filterBitSet;


    protected HalfByteGenomeSequenceGVCF(Map<Chromosome, byte[]>chromPositionMap, PositionList gvcfAnnotationsAndCalls) {
        this.chromPositionMap = chromPositionMap;
        this.gvcfAnnotationsAndCalls = gvcfAnnotationsAndCalls;
        chromPositionMap.entrySet().stream()
                .forEach(e -> chromLengthLookup.put(e.getKey(),e.getKey().getLength()));
        LongAdder genomeIndex=new LongAdder();
        chromosomes().stream().sorted()
                .forEach(chrom -> {
                    int length=chromLengthLookup.get(chrom);
                    wholeGenomeIndexMap.put(Range.closed(genomeIndex.longValue(),
                            genomeIndex.longValue()+length-1),chrom);
                    genomeIndex.add(length);}
                );
        genomeSize=genomeIndex.longValue();
        maskBitSet = new OpenBitSet(gvcfAnnotationsAndCalls.size());
        filterBitSet = new OpenBitSet(gvcfAnnotationsAndCalls.size());
    }
    @Override
    public Set<Chromosome> chromosomes() {
        return chromPositionMap.keySet();
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom) {
        return chromosomeSequence(chrom,1,chromLengthLookup.get(chrom));
    }

    @Override
    // Code assumes 1-based coordinates have been passed.  It will catch and return
    // null if the startSite is 0.  Otherwise, the user is on their own to ensure
    // input is 1-based.
    //TODO add in functionality using GVCF annotations
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int lastSite) {
        final int startSiteFinal = startSite;
        final int lastSiteFinal = lastSite;
        //zrm22 new for getting gvcf records.
        //Basically the loop now changes to loop through the positionList to check if the allele should be the ref or not
        gvcfAnnotationsAndCalls.chromosomeSiteCount(chrom);
        ArrayList<Position> listOfChrPositions = (ArrayList<Position>)gvcfAnnotationsAndCalls.stream()
                                                                        .filter(position->position.getChromosome().equals(chrom))
                                                                        .filter(position -> position.getPosition()>=startSiteFinal && position.getPosition()<=lastSiteFinal)
                                                                        .collect(Collectors.toList());
        //sort things to make sure we can iterate properly
        Collections.sort(listOfChrPositions);
        //loop through from startSite till the first position in the GVCF export these alleles directly from the ref
        int listOfChrPositionsCounter = 0;
        ArrayList<Byte> byteList = new ArrayList<>();
        for(int siteCounter = startSite; siteCounter <= lastSite && listOfChrPositionsCounter<listOfChrPositions.size(); siteCounter++) {
            if(listOfChrPositions.get(listOfChrPositionsCounter).getPosition()==siteCounter) {
                //check to see if the current position is a reference block or not
                SetMultimap<String,String> annos = listOfChrPositions.get(listOfChrPositionsCounter).getAnnotation().getAnnotationAsMap();
                if(annos.containsKey("END")) {
                    //this means we have a block
                    //Check the call and grab the corresponding allele values
                    String call = (String)annos.get("GT").toArray()[0];
                    boolean phased = true;
                    if(call.contains("/")) {
                        phased = false;
                    }
                    String[] callSplit = phased?call.split("|"):call.split("/");
                    int leftAllele = Integer.parseInt(callSplit[0]);
                    int rightAllele = Integer.parseInt(callSplit[1]);
                    int endPoint = Integer.parseInt((String)annos.get("END").toArray()[0])-1;
                    //Get the known Variants
                    int startSiteShifted = siteCounter - 1;  //shift over to zero base

                    if (startSite < 0) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: starting parameter is less than 1 for 1-based method");; // method needs 1-based coordinates.
                    byte[] packedBytes = chromPositionMap.get(chrom);
                    if (packedBytes == null) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); // chromosome not found
                    if (startSiteShifted > packedBytes.length*2 || endPoint > packedBytes.length*2 ) {
                        throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested sequence is out of range"); // requested sequence is out of range
                    }
                    String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                    String leftAlleleString = variants[leftAllele];
                    //TODO check to see if we should assume only Ref or<NON_REF> call for blocks
                    if(leftAllele == 0) {
                        //fill in with reference sequence
                        //pull the sequence from siteCounter till you get to endPoint

                        for (int i = startSiteShifted; i <= endPoint; i++) {
                            byteList.add((byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F)));
                        }

                    }
                    else {

                        for (int i = startSiteShifted; i <= endPoint; i++) {
                            byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                        }
                    }


                    //shift up the siteCounter to match the end
                    siteCounter=endPoint;
                }
                else {
                    //it is likely a snp
                    //check the call and grab the correct allele
                    //siteCounter will increment correctly

                    String call = (String)annos.get("GT").toArray()[0];
                    boolean phased = true;
                    if(call.contains("/")) {
                        phased = false;
                    }
                    String[] callSplit = phased?call.split("|"):call.split("/");
                    int leftAllele = Integer.parseInt(callSplit[0]);
                    int rightAllele = Integer.parseInt(callSplit[1]);
                    //Get the known Variants
                    String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                    //TODO handle hets correctly for indels and such
                    String leftAlleleString = variants[leftAllele];
                    for(int i = 0; i < leftAlleleString.length(); i++) {
                        byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                    }
                }


                listOfChrPositionsCounter++;
            }
            else {
                //grab the reference as we dont have a gvcf for the requested site
                //Should this be the breaking point of the sequence??
                //TODO should we mark these with Ns?
            }
        }
        byte[] fullBytes = new byte[byteList.size()];
        for(int i = 0; i < fullBytes.length; i++) {
            fullBytes[i] = byteList.get(i);
        }



//        startSite--;  //shift over to zero base
//        lastSite--;   //shift over to zero base
//        if (startSite < 0) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: starting parameter is less than 1 for 1-based method");; // method needs 1-based coordinates.
//        byte[] packedBytes = chromPositionMap.get(chrom);
//        if (packedBytes == null) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); // chromosome not found
//        if (startSite > packedBytes.length*2 || lastSite > packedBytes.length*2 ) {
//            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested sequence is out of range"); // requested sequence is out of range
//        }
//        byte[] fullBytes = new byte[lastSite - startSite + 1];
//        for (int i = startSite; i <= lastSite; i++) {
//            fullBytes[i - startSite] = (byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F));
//        }
        return fullBytes;
    }

    @Override
    //TODO add in functionality to handle GVCF annotations
    public byte[] genomeSequence(long startSite, long lastSite) {
        if(lastSite-startSite>Integer.MAX_VALUE) throw
                new IllegalArgumentException("Less than "+Integer.MAX_VALUE+" sites must be requested at a time");
        byte[] fullBytes=new byte[(int)(lastSite-startSite+1)];
        long currentSiteToGet=startSite;
        while(currentSiteToGet<lastSite) {
            Map.Entry<Range<Long>,Chromosome> rangeChromEntry=wholeGenomeIndexMap.getEntry(currentSiteToGet);
            int chrStart=(int)(currentSiteToGet-rangeChromEntry.getKey().lowerEndpoint());
            int chrLast=(int)Math.min(rangeChromEntry.getKey().upperEndpoint()-rangeChromEntry.getKey().lowerEndpoint(),lastSite-rangeChromEntry.getKey().lowerEndpoint());
            byte[] chromoSeq=chromosomeSequence(rangeChromEntry.getValue(), chrStart+1,chrLast+1);  //+1 for 0 based genome, 1 based chromosomes
            System.arraycopy(chromoSeq,0,fullBytes,(int)(currentSiteToGet-startSite),chromoSeq.length);
            currentSiteToGet+=chromoSeq.length;
        }
        return fullBytes;
    }

    @Override
    public int chromosomeSize(Chromosome chromosome) {
        return chromLengthLookup.get(chromosome);
    }

    @Override
    public long genomeSize() {
        return genomeSize;
    }

    @Override
    public int numberOfChromosomes() {
        return chromPositionMap.size();
    }

    @Override
    //TODO see if we need to fix this for GVCF annotations
    public Map<Long, Tuple<Chromosome, Integer>> fullRefCoordinateToChromCoordinate(ArrayList<Long> coordinates) {
        // Returns 0-based value from Chromosome array (values are stored as 0-based)
        Map<Long, Tuple<Chromosome, Integer>> mappedCoordinates = new ConcurrentHashMap<Long, Tuple<Chromosome, Integer>>();
        coordinates.stream().parallel().forEach(coordinate -> {
            Map.Entry<Range<Long>,Chromosome> rangeChromEntry=wholeGenomeIndexMap.getEntry(coordinate);
            Chromosome curChrom = rangeChromEntry.getValue();
            long chromCoordinate = coordinate - rangeChromEntry.getKey().lowerEndpoint();
            Tuple<Chromosome, Integer> chromWithCoordinate = new Tuple<>(curChrom, (int)chromCoordinate);
            mappedCoordinates.put(coordinate, chromWithCoordinate);
        });
        return mappedCoordinates;
    }

    public HashMap<Chromosome,ArrayList<ArrayList<Integer>>> getConsecutiveRegions() {
        HashMap<Chromosome,ArrayList<ArrayList<Integer>>> consecRegions = new HashMap<>();
        Set<Chromosome> chromosomeSet = chromosomes();
        //Loop through each chromosome add to a Map<RangeMap>
        HashMap<Chromosome,RangeSet<Integer>> rangeMaps = new HashMap<>();
        for(Chromosome chr : chromosomeSet) {
            rangeMaps.put(chr,TreeRangeSet.create());
            consecRegions.put(chr,new ArrayList<>());
        }

        for(int i = 0 ; i < gvcfAnnotationsAndCalls.size(); i++) {
            Position currentPosition = gvcfAnnotationsAndCalls.get(i);
            SetMultimap<String,String> annos = currentPosition.getAnnotation().getAnnotationAsMap();
            if(annos.containsKey("END")) {
                int endPoint = Integer.parseInt((String)annos.get("END").toArray()[0]);
                rangeMaps.get(currentPosition.getChromosome()).add(Range.closed(currentPosition.getPosition(),endPoint+1));
            }
            else {
                rangeMaps.get(currentPosition.getChromosome()).add(Range.closed(currentPosition.getPosition(),currentPosition.getPosition()+1));
            }

        }

        for(Chromosome chr : chromosomeSet) {
            Set<Range<Integer>> rangeSet = rangeMaps.get(chr).asRanges();
            for(Range<Integer> currentRange : rangeSet) {
                ArrayList<Integer> currentList = new ArrayList<Integer>();
                currentList.add(currentRange.lowerEndpoint());
                //We need to subtract 1 point so it will be [inclusive,inclusive] instead of [inclusive,exclusive)
                currentList.add(currentRange.upperEndpoint()-1);
                consecRegions.get(chr).add(currentList);
            }
        }

        return consecRegions;
    }
    //TODO move this to a different class
    public void writeFASTA(String fileName){
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
            HashMap<Chromosome,ArrayList<ArrayList<Integer>>> consecutiveRegions = getConsecutiveRegions();
            Set<Chromosome> chromosomes = chromosomes();
            for(Chromosome chr : chromosomes) {
                for(ArrayList<Integer> bounds : consecutiveRegions.get(chr)) {
                    writer.write(">Chr_"+chr.getChromosomeNumber()+"_StartSite_"+bounds.get(0)+"_EndSite_"+bounds.get(1));
                    writer.newLine();
                    writer.write(""+NucleotideAlignmentConstants.nucleotideBytetoString(chromosomeSequence(chr,bounds.get(0),bounds.get(1))));
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }

    }

    public BitSet getMaskBitSet() {
        return maskBitSet;
    }
    public void setMaskBitSet(BitSet newMaskBitSet) {
        maskBitSet = newMaskBitSet;
    }

    public BitSet getFilterBitSet() {
        return filterBitSet;
    }
    public void setFilterBitSet(BitSet newFilterBitSet) {
        filterBitSet = newFilterBitSet;
    }

    public void flipMaskBit(int index) {
        maskBitSet.fastFlip(index);
    }
    public void flipFilterBit(int index) {
        filterBitSet.fastFlip(index);
    }
}

class GVCFRecord {
    //Simple class to hold a gvcf record in a storage efficient manner
    //Need the following:
    //the call
    int[] call = new int[2];
    //phasing
    boolean phased = false;
    //the known variants
    ArrayList<String> knownVariants = new ArrayList<>();
    //the depth
    int depth = 0;
    //GQ score
    int gqScore = 0;
    //GeneralAnnotations for any other in INFO tag
    GeneralAnnotationStorage generalAnnotationStorage;


}

class GVCFPositionRecord {
    GeneralPosition ap= new GeneralPosition.Builder(new Chromosome("1"),1232)
            .maf(0.05f)
            .build();
    HeaderPositions hp;

    public GVCFPositionRecord(HeaderPositions hp) {
        this.hp = hp;

    }
    public Position parseGVCFRecords(String input) {
        //Figure out the tab positioning for the header columns
        int[] tabPos=new int[hp.NUM_HAPMAP_NON_TAXA_HEADERS+1];
        int tabIndex=0;
        int len=input.length();
        for (int i=0; (tabIndex<hp.NUM_HAPMAP_NON_TAXA_HEADERS+1)&&(i<len); i++) {
            if (input.charAt(i)=='\t') {
                tabPos[tabIndex++]=i;
            }
        }
        String chrName=input.substring(0, tabPos[hp.CHROMOSOME_INDEX]);
        Chromosome currChr=new Chromosome(new String(chrName));

        String snpID=null;
        if(hp.SNPID_INDEX>0) snpID=input.substring(tabPos[hp.SNPID_INDEX-1]+1, tabPos[hp.SNPID_INDEX]);

        //create the General position
        GeneralPosition.Builder apb=new GeneralPosition.Builder(currChr, Integer.parseInt(input.substring(tabPos[hp.POSITION_INDEX-1]+1, tabPos[hp.POSITION_INDEX])));
        if(snpID!=null && !snpID.equals(".")) {
            apb.snpName(snpID);
        }



        String refS=input.substring(tabPos[hp.REF_INDEX-1]+1, tabPos[hp.REF_INDEX]);
        String alt=input.substring(tabPos[hp.ALT_INDEX-1]+1, tabPos[hp.ALT_INDEX]);
        //create an String to hold the list of variants, we will have to parse this on export
        String allAlleles = refS+"/"+alt.replace(",","/");

        apb = apb.knownVariants(allAlleles);

        //loop through the INFO tag and save off any of those values into positions
        for(String annoS: Splitter.on(";").split(input.substring(tabPos[hp.INFO_INDEX-1]+1, tabPos[hp.INFO_INDEX]))) {
            apb.addAnno(annoS);
        }


        final int iGT=0; //genotype index
        int iAD=-1,iDP=-1,iGQ=-1, iPL=-1;  //alleleDepth, overall depth, genotypeQuality, phredGenotypeLikelihoods
        if(hp.FORMAT_INDEX>=0) {
            //Check to see if FORMAT tag is missing. Only applicable for single taxa files
            if(tabPos[hp.FORMAT_INDEX]==0) {
                throw new IllegalStateException("Error Processing VCF: Missing FORMAT tag.");
            }
            String unsplitInput = input.substring(tabPos[hp.FORMAT_INDEX-1]+1, tabPos[hp.FORMAT_INDEX]);
            if(unsplitInput.length()==0|| !unsplitInput.startsWith("GT")) {
                //Check to see it has the GT field
                if(unsplitInput.contains("GT")) {
                    throw new IllegalStateException("Error Processing VCF Block: GT field is not in first position of FORMAT.");
                }
                //If GT isnt in, we assume that it is missing FORMAT
                else {
                    throw new IllegalStateException("Error Processing VCF Block: Missing FORMAT tag.");
                }
            }
            String[] formatS = unsplitInput.split(":");

            iAD=firstEqualIndex(formatS,"AD");
            iDP=firstEqualIndex(formatS,"DP");
            iGQ=firstEqualIndex(formatS,"GQ");
        }



        //after info is recorded we need to record the call section
        String taxaAllG = input.substring(tabPos[hp.NUM_HAPMAP_NON_TAXA_HEADERS-1]+1);
        int f = 0;
        for(String fieldS: Splitter.on(":").split(taxaAllG)) {
            if (f == iGT) {
                apb.addAnno("GT",fieldS);
            }
            if(f == iAD) {
                apb.addAnno("AD",fieldS);
            }
            if(f == iDP) {
                apb.addAnno("DP",fieldS);
            }
            if(f == iGQ) {
                apb.addAnno("GQ",fieldS);
            }

            f++;
        }



        return apb.build();
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i=0; i<sa.length; i++) {
            if(sa[i].equals(match)) return i;
        }
        return -1;
    }
}