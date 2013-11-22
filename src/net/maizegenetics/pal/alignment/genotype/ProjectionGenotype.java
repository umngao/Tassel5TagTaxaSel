package net.maizegenetics.pal.alignment.genotype;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.dna.snp.AlignmentUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.DonorHaplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.NavigableSet;

/**
 * Projection genotype use defined haplotypes and breakpoints that point to high density
 * genotypes (baseAlignment).  These are used to efficiently store and connect low density maps with imputed high density genotypes.
 * <p></p>
 * The alignment built by this builder is a CoreAlignment with a ProjectionGenotype.  The taxa indices come from the
 * projection alignment file, while the site indices are the same as the base alignment.
 * TODO this implement a projection interface with the getDonorHaplotypes method
 *
 * @author Ed Buckler
 */
public class ProjectionGenotype extends AbstractGenotype {
    private final Alignment myBaseAlignment;  //high density marker alignment that is being projected. It was suggested that this
    //just have a pointer to a genotype, which would work, excepting for saving the file, when the base taxa names are needed.
    private ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints;

    private enum BaseMode {General, Site, Taxa};
    private BaseMode currMode=BaseMode.Taxa;
    private ArrayList<RangeMap<Integer,DonorSiteHaps>> breakMaps;
    private DonorSiteHaps[] currentDSH;
    private byte[] donorForCachedSite;
    private byte[] projForCachedTaxon;
    private int cachedSite=-1;
    private int cachedTaxon=-1;
    int[] primDSH; //startSite,endSite,parent1,parent2 array for the

    public ProjectionGenotype(Alignment hdAlign, ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints) {
        super(allBreakPoints.size(), hdAlign.getSiteCount(), false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myBaseAlignment = hdAlign;
        this.allBreakPoints=allBreakPoints;
        breakMaps=new ArrayList<>(getTaxaCount());
        for (NavigableSet<DonorHaplotypes> allBreakPoint : allBreakPoints) {
            RangeMap<Integer,DonorSiteHaps> tRM=TreeRangeMap.create();
            for (DonorHaplotypes dh : allBreakPoint) {
                int[] siteRange=siteRangeForDonor(dh);
                DonorSiteHaps dsh=new DonorSiteHaps(siteRange[0],siteRange[1],dh.getParent1index(),dh.getParent2index());
                tRM.put(Range.closed(siteRange[0],siteRange[1]),dsh);
                //TODO consider putting in blank range maps
            }
            breakMaps.add(tRM);
        }
        currentDSH=new DonorSiteHaps[getTaxaCount()];
        primDSH=new int[myTaxaCount*4];
        Arrays.fill(primDSH,Integer.MIN_VALUE);
    }

    public NavigableSet<DonorHaplotypes> getDonorHaplotypes(int taxon) {
        return allBreakPoints.get(taxon);
    }

    private int[] siteRangeForDonor(DonorHaplotypes dh) {
        int start=myBaseAlignment.getSiteOfPhysicalPosition(dh.getStartPosition(),dh.getChromosome());
        if(start<0) start=-(start+1);
        int end=myBaseAlignment.getSiteOfPhysicalPosition(dh.getEndPosition(),dh.getChromosome());
        if(end<0) end=-(end+1);
        return new int[]{start,end};
    }

    @Override
    public byte getBase(int taxon, int site) {
        if(currMode==BaseMode.Site) {return getBaseSite(taxon, site);}
//        if(currMode==BaseMode.Taxa) {return getBaseTaxon(taxon, site);}
        return getBaseGeneral(taxon, site);
    }

    /**
     * Returns the high density base alignment of the projection genotype.
     * @return base Alignment
     */
    public Alignment getBaseAlignment() {
        return myBaseAlignment;
    }

    private byte getBaseGeneral(int taxon, int site) {
        if((currentDSH[taxon]==null)||(!currentDSH[taxon].containsSite(site))) {
            currentDSH[taxon]=breakMaps.get(taxon).get(site);
            if(currentDSH[taxon]==null) return Alignment.UNKNOWN_DIPLOID_ALLELE;
            //TODO consider null
        }
        byte p1=myBaseAlignment.getBase(currentDSH[taxon].getParent1index(),site);
        byte p2=myBaseAlignment.getBase(currentDSH[taxon].getParent2index(),site);
        return AlignmentUtils.getUnphasedDiploidValueNoHets(p1, p2);
    }

    //Currently this is no faster than general getBase.  it should be possible to make this faster.
    private byte getBaseTaxon(int taxon, int site) {
        if(taxon!=cachedTaxon) {
            projForCachedTaxon=new byte[mySiteCount];
            Arrays.fill(projForCachedTaxon,Alignment.RARE_DIPLOID_ALLELE);
            cachedTaxon=taxon;
        }
        byte result=projForCachedTaxon[site];
        if(result==Alignment.RARE_DIPLOID_ALLELE) {
            DonorSiteHaps currentDSH=breakMaps.get(taxon).get(site);
            byte[] r=myBaseAlignment.getBaseRange(currentDSH.getParent1index(),currentDSH.getStartSite(),currentDSH.getEndSite()+1);
            System.arraycopy(r,0,projForCachedTaxon,currentDSH.getStartSite(),r.length);
            result=projForCachedTaxon[site];
        }
        return result;
    }

    private byte getBaseSite(int taxon, int site) {
        //test transpose problems
        if(site!=cachedSite) {
            donorForCachedSite=myBaseAlignment.getGenotypeMatrix().getGenotypeForAllTaxa(site);
            cachedSite=site;
        }
        int primPos=taxon<<2;
        if((site<primDSH[primPos++])||(site>primDSH[primPos++])) {
            DonorSiteHaps currentDSH=breakMaps.get(taxon).get(site);
            primPos=taxon<<2;
            primDSH[primPos++]=currentDSH.getStartSite();
            primDSH[primPos++]=currentDSH.getEndSite();
            primDSH[primPos++]=currentDSH.getParent1index();
            primDSH[primPos]=currentDSH.getParent2index();
            primPos=(taxon<<2)+2;
            //TODO consider null
        }
        //       if(primDSH[primPos]==primDSH[primPos+1]) return donorForCachedSite[primDSH[primPos]];
        return AlignmentUtils.getUnphasedDiploidValueNoHets(donorForCachedSite[primDSH[primPos]], donorForCachedSite[primDSH[primPos+1]]);
    }



    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }


    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBaseAlignment.getGenotypeMatrix().transposeData(siteInnerLoop);
        if(siteInnerLoop) {currMode=BaseMode.Site;}
        else {currMode=BaseMode.General;}

    }


   private class DonorSiteHaps {
        private final int startSite;
        private final int endSite;
        private final int parent1index;
        private final int parent2index;

       private DonorSiteHaps(int startSite, int endSite, int parent1index, int parent2index) {
           this.startSite=startSite;
           this.endSite=endSite;
           this.parent1index=parent1index;
           this.parent2index=parent2index;
       }

       private int getStartSite() {
           return startSite;
       }

       private int getEndSite() {
           return endSite;
       }

       private int getParent1index() {
           return parent1index;
       }

       private int getParent2index() {
           return parent2index;
       }

       private boolean containsSite(int site) {
           if((site<startSite)||(site>endSite)) return false;
           return true;
       }
   }

}
