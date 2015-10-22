package net.maizegenetics.analysis.rna;

import com.google.common.collect.*;
import com.google.common.primitives.Ints;
import net.maizegenetics.dna.tag.Tag;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Find Match by Kmers using a simple hashmap of fixed length kmers to determine most similar sequence.
 *
 * @author Ed Buckler
 * @author Karl Kremling
 *
 */
public class FindMatchByKmers {
    public enum MatchType{COMPLETE_SEQ, FIRST_KMER, MODAL_KMER};
    private final List<Tag> tags;
    private final Map<String,int[]> kmerMap;
    private final Map<String,Integer> completeSeqMap;
    private final int kmerLength;
    private final int maxKmerCopies;
    private final MatchType matchType;
    private final int[] emptyInts=new int[0];

    /**
     * Create kmer
     * @param tagSet unique set sequences to be indexed with kmers
     * @param kmerLength length of all the kmers
     * @param maxKmerCopies retain only those kmers appearing less then this value
     */
    public FindMatchByKmers(Set<Tag> tagSet, MatchType matchType, int kmerLength, int maxKmerCopies) { //find a contig which is the consensus source of the kmers from a read
        this.kmerLength=kmerLength; // length of kmers for matching
        this.maxKmerCopies=maxKmerCopies; // kmer appearance frequency cutoff. Kmers are later discarded if they appear in the contigs more than maxKmerCopies times
        this.matchType=matchType;
        tags=new ArrayList<>(tagSet); //holds all the sequences or contigs in an array
        long bpLength=0;
        Multimap<String,Integer> completeKmerMap = HashMultimap.create(4_000_000,2); //? hashmap pointing from each kmer to its Contig of origin
        for (int ti = 0; ti < tags.size(); ti++) { //iterates through the ArrayList of tags
            String seq=tags.get(ti).sequence(); //pulls the sequence of contig/tag ti
            bpLength+=seq.length(); //increase the bpLength counter by ti's length
            for (int i = 0; i <= seq.length() - kmerLength; i++) { //slides through the contig to pull out each kmer
                String sub = seq.substring(i, i + kmerLength);
                completeKmerMap.put(sub, ti); //appends a key value pair to the KmerMap in which the k=kmer and the v=contig of origin
            }
        }
        System.out.printf("Total bp %,d Kmer Map Size: %d %n", bpLength,completeKmerMap.size());
        //purge high copy kmers from the multimap
        //consider making it immutable
        kmerMap=new HashMap<>(completeKmerMap.size()*3/2); //? why make it 1.5x larger than the complete kmer multimap if we know we will have kmers <= CompleteKmerMap.size?
        for (String kmer : completeKmerMap.keySet()) {
            Collection<Integer> tagIndices=completeKmerMap.get(kmer);
            if(tagIndices.size()< this.maxKmerCopies) { //only record key-value pairs for kmers appearing less than maxKmerCopies times
                kmerMap.put(kmer, Ints.toArray(tagIndices));
            }
        }
        //Allow direct searching for the complete sequence
        completeSeqMap = new HashMap<>(tags.size()*3/2);
        for (int ti = 0; ti < tags.size(); ti++) {
            completeSeqMap.put(tags.get(ti).sequence(),ti);
        }
        System.out.println(kmerMap.size());

    }

    public int totalNumberOfTags() {
        return tags.size();
    }

    public List<Tag> getTags() {
        return Collections.unmodifiableList(tags);
    }

    /**
     * Return the first unique kmer match found starting at the 5' end of the sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    public Optional<Tag> getMatchTag(String seq) {
        OptionalInt tagIndex=getMatchIndex(seq);
        return tagIndex.isPresent()?Optional.of(tags.get(tagIndex.getAsInt())):Optional.empty();
    }

    /**
     * Return the first unique kmer match found starting at the 5' end of the sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    public OptionalInt getMatchIndex(String seq) {
        switch(matchType) {
            case COMPLETE_SEQ:return lookForCompleteMatch(seq);
            case FIRST_KMER:return getFirstUniqueMatchIndex(seq);
            case MODAL_KMER:return getMostCommonMatchIndex(seq);
            default:
                System.err.println("Bad Match Type");
                return OptionalInt.empty();
        }
    }

    /**
     * Return the first unique kmer match found starting at the 5' end of the sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    private OptionalInt getFirstUniqueMatchIndex(String seq) {
        return IntStream.range(0,seq.length() - kmerLength)
                .mapToObj(i -> kmerMap.getOrDefault(seq.substring(i, i + kmerLength), emptyInts)) //search the kmer map
                .filter(hits -> hits.length == 1)  //only return kmer hits with a single match
                .mapToInt(hits -> hits[0])
                .findFirst();
    }

    /**
     * Returns the sequence match with the most kmer hits across the entire length of the query sequence.
     * @param seq query sequence - usually the reads
     * @return matching sequence
     */
    private OptionalInt getMostCommonMatchIndex(String seq) {
        int[] allHits=IntStream.range(0, seq.length() - kmerLength).filter(i -> i%3==0) //stream formatting of a loop to only consider kmers starting at evey 3rd position
                .flatMap(i -> IntStream.of(kmerMap.getOrDefault(seq.substring(i, i + kmerLength), emptyInts))) // load kmer from pos i through kmerLength
                .toArray();
        if(allHits.length==0) return OptionalInt.empty(); //if there are no hits for any of the kmers in a read, return an emtpy hit entry
        return OptionalInt.of(mode(allHits)); // returns the most frequently hit contig for all the kmers in the read
    }

    private OptionalInt lookForCompleteMatch(String seq) {
        Integer tagIndex=completeSeqMap.get(seq);
        if (tagIndex!=null) {return OptionalInt.of(tagIndex);}
        else {return OptionalInt.empty();}
    }

    private static int mode(int[] array) {
        HashMap<Integer, Integer> cntMap = new HashMap<>();
        int maxCnt = 1, mode = array[0];
        for (int i = 0; i < array.length; i++) {
            int count = cntMap.compute(array[i], (k, cnt) -> cnt == null ? 1 : cnt + 1);
            if (count > maxCnt) {
                maxCnt = count;
                mode = array[i];
            }
        }
        return mode;
    }



}
