package net.maizegenetics.dna.snp.io;

import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 * Create an alignment based on HapMap format file (either .txt or compressed).  Alleles are set as global major and
 * global minor.
 * e.g. code <p></p>
 * {@code
 * Alignment a=BuilderFromHapMap.getBuilder(infileName).build();
 * }
 * <p></p>
 * TODO:  Add filtering while reading, provide an option to define the alleles as reference and alternate
 *
 * @author Ed Buckler
 */
public class BuilderFromHapMap {

    private static final Logger myLogger=Logger.getLogger(BuilderFromHapMap.class);
    private static final Pattern WHITESPACE_PATTERN=Pattern.compile("\\s");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS=11;
    private final String infile;
    private int[] taxaRedirect;
    private boolean sortAlphabetically=false;

    private BuilderFromHapMap(String infile) {
        this.infile=infile;
    }

    public static BuilderFromHapMap getBuilder(String infile) {
        return new BuilderFromHapMap(infile);
    }

    //TODO provide options on caching to use, read only some sites, etc.
    public GenotypeTable build() {
        long time=System.nanoTime();
        GenotypeTable result=null;
        try {
            int numThreads=Runtime.getRuntime().availableProcessors();
            ExecutorService pool=Executors.newFixedThreadPool(numThreads);
            System.out.print("Reading :"+infile+" ");
            BufferedReader r=Utils.getBufferedReader(infile, -1);
            Map<String,SetMultimap<String,String>> sampAnnoBuild=new TreeMap<>();
            String currLine;
            while (((currLine=r.readLine())!=null)&&(currLine.startsWith("##"))) {
                String[] cat=currLine.split("=",2);
                if(cat.length<2) continue;
                if(cat[0].startsWith("##SAMPLE")) {
                    SetMultimap<String, String> mapOfAnno=TaxaListIOUtils.parseVCFHeadersIntoMap(cat[1]);
                    String taxaID=mapOfAnno.get("ID").iterator().next();
                    if(taxaID==null) break;
                    sampAnnoBuild.put(taxaID,mapOfAnno);
                }
            }
            TaxaList taxaList=processTaxa(currLine,sampAnnoBuild);
            int linesAtTime=1<<12;  //this is a critical lines with 20% or more swings.  Needs to be optimized with transposing
            //  int linesAtTime=1<<8;  //better for with lots of taxa.
            ArrayList<String> commentLines=new ArrayList<>();
            ArrayList<String> txtLines=new ArrayList<>(linesAtTime);
            ArrayList<ProcessHapMapBlock> pbs=new ArrayList<>();
            int lines=0;
            while ((currLine=r.readLine())!=null) {
                txtLines.add(currLine);
                lines++;
                if (lines%linesAtTime==0) {
                    ProcessHapMapBlock pb=ProcessHapMapBlock.getInstance(taxaList.numberOfTaxa(), txtLines);
                    pbs.add(pb);
                    //     pb.run();
                    pool.execute(pb);
                    txtLines=new ArrayList<>(linesAtTime);
                    System.out.print(".");
                }
            }
            r.close();
            if (txtLines.size()>0) {
                ProcessHapMapBlock pb=ProcessHapMapBlock.getInstance(taxaList.numberOfTaxa(), txtLines);
                pbs.add(pb);
                pool.execute(pb);
            }
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BuilderFromHapMap: processing threads timed out.");
            }
            int currentSite=0;
            PositionListBuilder posBuild=new PositionListBuilder();
            taxaList=sortTaxaListIfNeeded(taxaList);
            GenotypeCallTableBuilder gb=GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.numberOfTaxa(), lines);
            for (ProcessHapMapBlock pb : pbs) {
                posBuild.addAll(pb.getBlkPosList());
                SuperByteMatrix bgTS=pb.getGenoTS();
                for (int t=0; t<bgTS.getNumRows(); t++) {
                    int currentTaxa = taxaRedirect[t];
                    for (int s=0; s<bgTS.getNumColumns(); s++) {
                        gb.setBase(currentTaxa, currentSite+s, bgTS.get(t, s));
                    }
                }
                currentSite+=pb.getSiteNumber();
            }
            if (posBuild.validateOrdering()==false) {
                throw new IllegalStateException("BuilderFromHapMap: Ordering incorrect HapMap must be ordered by position");
            }
            GenotypeCallTable g=gb.build();
            result=GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList);
            System.out.println(" finished");
        } catch (IOException|InterruptedException e) {
            e.printStackTrace();
        }
        long totalTime=System.nanoTime()-time;
       // System.out.printf("BuilderFromHapMap data timing %gs %n", totalTime/1e9);
        return result;
    }

    /*
   Set the builder so that when built it will sort the taxa
    */
    public BuilderFromHapMap sortTaxa() {
        sortAlphabetically=true;
        return this;
    }

    private TaxaList sortTaxaListIfNeeded(TaxaList origTL) {
        TaxaList resultTL;
        if(sortAlphabetically){
            resultTL=new TaxaListBuilder().addAll(origTL).sortTaxaAlphabetically().build();
        } else {
            resultTL=origTL;
        }
        taxaRedirect=new int[origTL.numberOfTaxa()];
        for (int i=0; i<taxaRedirect.length; i++) {
            taxaRedirect[i]=resultTL.indexOf(origTL.get(i));
        }
        return resultTL;
    }

    private static TaxaList processTaxa(String readLn, Map<String,SetMultimap<String,String>> taxaAnnotation) {
        String[] header=WHITESPACE_PATTERN.split(readLn);
        int numTaxa=header.length-NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i=0; i<numTaxa; i++) {
            String taxonID=header[i+NUM_HAPMAP_NON_TAXA_HEADERS];
            Taxon.Builder at=new Taxon.Builder(taxonID);
            SetMultimap<String,String> taMap=taxaAnnotation.get(taxonID);
            if(taMap!=null) {
                for (Map.Entry<String,String> en : taMap.entries()) {
                    if(en.getKey().equals("ID")) continue; //skip the IDs as these became the name
                    String s=en.getValue().replace("\"","");
                    at.addAnno(en.getKey(),s);
                }
            }
            tlb.add(at.build());
        }
        return tlb.build();
    }
}

class ProcessHapMapBlock implements Runnable {

    private static final Pattern WHITESPACE_PATTERN=Pattern.compile("\\s");
    private static final Pattern SLASH_PATTERN=Pattern.compile("/");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS=11;
    private static final int GENOIDX=NUM_HAPMAP_NON_TAXA_HEADERS;
    private static final int SNPID_INDEX=0;
    private static final int VARIANT_INDEX=1;
    private static final int CHROMOSOME_INDEX=2;
    private static final int POSITION_INDEX=3;
    private final int taxaN;
    private final int siteN;
    private ArrayList<String> txtL;
    private SuperByteMatrix gTS;
    private final ArrayList<Position> blkPosList;
 //   private final byte[] convert;
    private final boolean isOneLetter; //true e.g. A,R, false=AA,CT

    private ProcessHapMapBlock(int taxaN, ArrayList<String> txtL) {
        this.taxaN=taxaN;
        this.siteN=txtL.size();
        this.txtL=txtL;
        blkPosList=new ArrayList<>(siteN);
        String[] tokens=WHITESPACE_PATTERN.split(txtL.get(0), NUM_HAPMAP_NON_TAXA_HEADERS+1);
        double avg=(double) (tokens[GENOIDX].length()+1)/(double) taxaN;
        if ((avg>1.99)&&(avg<2.01)) {
            isOneLetter=true;
        } else if ((avg>2.99)&&(avg<3.01)) {
            isOneLetter=false;
        } else {
            throw new IllegalStateException("ProcessHapMapBlock: Genotype coded wrong use 1 or 2 letters per genotype");
        }
    }

    public static ProcessHapMapBlock getInstance(int taxaN, ArrayList<String> txtL) {
        return new ProcessHapMapBlock(taxaN, txtL);
    }

    //@Override
    public void run() {
        Map<String, Chromosome> chromosomeLookup=new HashMap<>();
        gTS = SuperByteMatrixBuilder.getInstance(taxaN, siteN);
        for (int s=0; s<siteN; s++) {
            String input = txtL.get(s);
            try {
                int[] tabPos = new int[NUM_HAPMAP_NON_TAXA_HEADERS];
                int tabIndex = 0;
                int len = input.length();
                for (int i = 0; (tabIndex < NUM_HAPMAP_NON_TAXA_HEADERS) && (i < len); i++) {
                    if (input.charAt(i) == '\t') {
                        tabPos[tabIndex++] = i;
                    }
                }
                String chrName = input.substring(tabPos[CHROMOSOME_INDEX - 1] + 1, tabPos[CHROMOSOME_INDEX]);
                Chromosome currChr = chromosomeLookup.get(chrName);
                if (currChr == null) {
                    currChr = new Chromosome(new String(chrName));
                    chromosomeLookup.put(chrName, currChr);
                }
                String variants = input.substring(tabPos[VARIANT_INDEX - 1] + 1, tabPos[VARIANT_INDEX]);
                GeneralPosition.Builder apb = new GeneralPosition.Builder(currChr, Integer.parseInt(input.substring(tabPos[POSITION_INDEX - 1] + 1, tabPos[POSITION_INDEX])))
                        .snpName(input.substring(0, tabPos[SNPID_INDEX]))
                        .knownVariants(variants) //TODO                    strand, variants,
                        ;
                try {
                    byte glbMajor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(0));
                    apb.allele(WHICH_ALLELE.GlobalMajor, glbMajor);
                    if (variants.length() == 3) {
                        byte glbMinor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(2));
                        apb.allele(WHICH_ALLELE.GlobalMinor, glbMinor);
                    }
                } catch (IllegalArgumentException e) {
                    //for the indels that cannot be converted correctly now
                    // System.out.println("Error Parsing this variant"+Arrays.toString(variants));
                }
                blkPosList.add(apb.build());
                int offset = tabPos[NUM_HAPMAP_NON_TAXA_HEADERS - 1] + 1;
                if (isOneLetter) {
                    for (int i = offset; i < len; i += 2) {
                        gTS.set((i - offset) / 2, s, NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i)));
                    }
                } else {
                    for (int i = offset; i < len; i += 3) {
                        //there is a phasing conflict with the existing import approach
                        gTS.set((i - offset) / 3, s, GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i + 1)),
                                NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i))));
                    }
                }
                swapSitesIfOutOfOrder(s);
            } catch(Exception e) {
                System.err.println("Error parsing this row "+input);
                throw e;
            }
        }
        txtL=null;
    }

    //Swap adjacent misordered sites, often caused by two sites at the same positions with a different name order
    private void swapSitesIfOutOfOrder(int site) {
        if(site<1) return;
        if(blkPosList.get(site-1).compareTo(blkPosList.get(site))>0) {
            //swap
            Position tempP=blkPosList.get(site-1);
            System.out.println();
            System.out.print("Swapping:" + tempP.toString() + " <-> " + blkPosList.get(site).toString());
            blkPosList.set(site - 1, blkPosList.get(site));
            blkPosList.set(site,tempP);
            for (int t = 0; t < gTS.getNumRows(); t++) {
                byte tempG=gTS.get(t, site-1);
                gTS.set(t, site-1,gTS.get(t,site));
                gTS.set(t, site,tempG);
            }
        }

    }

    int getSiteNumber() {
        return siteN;
    }

    SuperByteMatrix getGenoTS() {
        return gTS;
    }

    ArrayList<Position> getBlkPosList() {
        return blkPosList;
    }
}
