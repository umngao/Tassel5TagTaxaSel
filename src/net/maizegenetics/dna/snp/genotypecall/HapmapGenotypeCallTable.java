package net.maizegenetics.dna.snp.genotypecall;

import gnu.trove.list.array.TLongArrayList;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.concurrent.ExecutionException;
import java.util.regex.Pattern;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

import org.apache.commons.lang.NotImplementedException;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;

public class HapmapGenotypeCallTable extends AbstractGenotypeCallTable {
    private static Pattern tab = Pattern.compile("\t");
    private FileChannel myGenotypeFileChannel; 
    private long[] linePointers;
    private static final Charset myCharset = Charset.forName("UTF-8");
    private static final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
    
    private LoadingCache<Integer, String[]> mySiteCache = CacheBuilder.newBuilder().maximumSize(100).build(new CacheLoader<Integer, String[]>(){

        @Override
        public String[] load(Integer arg) throws Exception {
            return readLine(arg);
        }
        
    } );
    
    public HapmapGenotypeCallTable(int numTaxa, long[] linePointers, String HapmapFilename) {
        super(numTaxa, linePointers.length - 1, false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        try {
            myGenotypeFileChannel = FileChannel.open(Paths.get(HapmapFilename), StandardOpenOption.READ);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(String.format("Unable to open %s for random access", HapmapFilename), e);
        } catch (IOException e) {
            throw new RuntimeException(String.format("Unable to open %s for random access", HapmapFilename), e);
        }
        this.linePointers = linePointers;
    }
    
    private String[] readLine(int site) throws IOException {
        long startLine = linePointers[site];
        long endLine = linePointers[site + 1];
        int size = (int) (endLine - startLine);
        ByteBuffer linebytes = ByteBuffer.allocate(size);
        int nread = myGenotypeFileChannel.read(linebytes, startLine);
        String input = new String(linebytes.array(), myCharset);
        return tab.split(input);
    }
    
    @Override
    public byte genotype(int taxon, int site) {
        String[] myLine;
        try {
            myLine = mySiteCache.get(site);
        } catch (ExecutionException e) {
            myLine = null;
        }
        
        if (myLine == null) return NN;
        return NucleotideAlignmentConstants.getNucleotideDiploidByte(myLine[taxon + 11]);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        String[] myLine;
        try {
            myLine = mySiteCache.get(site);
        } catch (ExecutionException e) {
            myLine = null;
        }
        
        if (myLine == null) return "N";
        return myLine[taxon + 11];
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new NotImplementedException("method transposeData not implemented in HapmapGenotypeCallTable.");
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        String[] info;
        try {
            info = readLine(site);
        } catch (IOException e) {
            throw new RuntimeException(String.format("Error getting hapmap record for site %d: ", site), e);
        }
        
        int ntaxa = numberOfTaxa();
        byte[] geno  = new byte[ntaxa];
        for (int t = 0; t < ntaxa; t++) geno[t] = NucleotideAlignmentConstants.getNucleotideDiploidByte(info[t + 11]);
        return geno;
    }
    
    public static GenotypeTable getInstanceforRandomAccessHapmapFile(String hapmapFile) {
        final byte eol = (byte) '\n';
        final byte cr = (byte) '\r';

        final Pattern tab = Pattern.compile("\t");
        TLongArrayList linePointers = new TLongArrayList();
        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        PositionListBuilder posBuilder = new PositionListBuilder();
        
        long start = System.nanoTime();
        int bufferSize = 128000;
        byte[] buffer = new byte[bufferSize];
        ByteBuffer mybb = ByteBuffer.wrap(buffer);
        long currentFilePos = 0;
        long filesize;
        long[] linePtr;
        try (FileChannel fc = (FileChannel.open(Paths.get(hapmapFile), StandardOpenOption.READ))) {
            filesize = fc.size();
            int nread;
            boolean endCRLF = false;
            do {
                mybb.clear();
                nread = fc.read(mybb);
                for (int i = 0; i < nread; i++) {
                    if (endCRLF) { //checks to see if buffer begins with eol
                        if (buffer[i] == eol) i++;
                        long pos = currentFilePos + i;
                        if (pos < filesize) linePointers.add(currentFilePos + i);
                        endCRLF = false;
                    }
                    if (buffer[i] == eol || buffer[i] == cr) {
                        i++;
                        if (i < nread && buffer[i] == eol) i++;
                        if (i < nread) {
                            long pos = currentFilePos + i;
                            if (pos < filesize) linePointers.add(currentFilePos + i);
                        } else endCRLF = true;
                    }
                }
                currentFilePos += nread;
            } while (nread == bufferSize);
            linePointers.add(filesize);
            
            System.out.printf("Read positions in %d ms.\n", (System.nanoTime() - start)/1000000);
            
            //read taxa
            start = System.nanoTime();
            Charset myCharset = Charset.forName("UTF-8");  
            linePtr = linePointers.toArray();
            int nlines = linePtr.length - 1;
            int size = (int) linePtr[0];
            ByteBuffer linebytes = ByteBuffer.allocate(size);
            nread = fc.read(linebytes, 0);
            String input = new String(linebytes.array(), myCharset);
            input = input.substring(0, input.length() - 1);
            if (input.endsWith("\r")) input = input.substring(0, input.length() - 1);
            String[] data = tab.split(input);
            String[] taxanames = Arrays.copyOfRange(data, 11, data.length);
            for (String tname : taxanames) 
                taxaBuilder.add(new Taxon(tname));
            
            System.out.printf("Created taxaList in %d ms.\n", (System.nanoTime() - start)/1000000);
            
            //read positions
            start = System.nanoTime();
            for (int i = 0; i < nlines; i++) {
                long startLine = linePtr[i];
                long endLine = linePtr[i +1];
                size = (int) (endLine - startLine);
                linebytes = ByteBuffer.allocate(size);
                nread = fc.read(linebytes, startLine);
                input = new String(linebytes.array(), myCharset);
                data = tab.split(input, 5);
                posBuilder.add(new GeneralPosition.Builder(new Chromosome(data[2]), Integer.parseInt(data[3])).snpName(data[0]).build());
            }
            System.out.printf("Created PositionList in %d ms.\n", (System.nanoTime() - start)/1000000);
        } catch (IOException ioe) {
            throw new RuntimeException("Error reading hapmap file: ", ioe);
        }
        
        System.out.println("finished reading in the positions");
        System.out.printf("%d positions, elapsed time = %d.\n", linePointers.size() - 1, (System.nanoTime() - start)/1000000);
        TaxaList myTaxa = taxaBuilder.build();
        return GenotypeTableBuilder.getInstance(new HapmapGenotypeCallTable(myTaxa.size(), linePtr,hapmapFile), posBuilder.build(), myTaxa);
    }
    
}
