package net.maizegenetics.analysis.imputation;

import com.google.common.base.Joiner;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.MinMaxPriorityQueue;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.analysis.popgen.LDResult;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.SimpleTableReport;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Need to fill in this
 *
 * @author Ed Buckler
 */
public class ImputationAccuracyPlugin extends AbstractPlugin {


    private static final Logger myLogger = Logger.getLogger(ImputationAccuracyPlugin.class);

    public ImputationAccuracyPlugin() {
        super(null, false);
    }

    public ImputationAccuracyPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }


    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.size() != 3) {
            throw new IllegalArgumentException("LDKNNiImputationPlugin: preProcessParameters: Please select three Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        // Load in the genotype table
        GenotypeTable origGenoTable = (GenotypeTable)input.getDataOfType(GenotypeTable.class).get(0).getData();
        GenotypeTable maskGenoTable = (GenotypeTable)input.getDataOfType(GenotypeTable.class).get(1).getData();
        GenotypeTable impGenoTable = (GenotypeTable)input.getDataOfType(GenotypeTable.class).get(2).getData();


        int[][] cnts=new int[5][5];
        for (int site = 0; site < origGenoTable.numberOfSites(); site++) {
            byte majorAllele=origGenoTable.majorAllele(site);
            byte minorAllele=origGenoTable.minorAllele(site);
            Map<Byte,Integer> genotypeToIndexMap=new HashMap<>();
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele,majorAllele),0);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(majorAllele,minorAllele),1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele,majorAllele),1);
            genotypeToIndexMap.put(GenotypeTableUtils.getDiploidValue(minorAllele,minorAllele),2);
            genotypeToIndexMap.put(GenotypeTable.UNKNOWN_DIPLOID_ALLELE,3);
            for (int taxonIdx = 0; taxonIdx < origGenoTable.numberOfTaxa(); taxonIdx++) {
                byte origGeno=origGenoTable.genotype(taxonIdx,site);
                if(origGeno==GenotypeTable.UNKNOWN_DIPLOID_ALLELE) continue;  //skip if unknown in the original data
                if(maskGenoTable.genotype(taxonIdx,site)!=GenotypeTable.UNKNOWN_DIPLOID_ALLELE) continue;  //skip if not missing in the masked data
                int originalIndex= genotypeToIndexMap.getOrDefault(origGeno,4);
                int impIndex= genotypeToIndexMap.getOrDefault(impGenoTable.genotype(taxonIdx, site), 4);
                if(impIndex==4) System.out.println(site+":"+impGenoTable.taxa().get(taxonIdx).toString()+":"+impGenoTable.genotype(taxonIdx, site) + NucleotideAlignmentConstants.getNucleotideIUPAC(impGenoTable.genotype(taxonIdx, site)));
                cnts[originalIndex][impIndex]++;
            }
        }

        //System.out.println(Arrays.deepToString(cnts));
        String[] headers={"AA","Aa","aa","N","Other"};
        int errors=0, correct=0;
        System.out.println(Joiner.on('\t').join(headers));
        for (int i = 0; i < cnts.length; i++) {
            System.out.print(headers[i]+"\t");
            System.out.println(Ints.join("\t", cnts[i]));
            correct+=cnts[i][i];
            errors+=cnts[i][0]+cnts[i][1]+cnts[i][2]-cnts[i][i];
        }
        double errorRate=(double)errors/(double)(correct+errors);
        System.out.printf("Correct:%d Errors:%d ErrorRate:%g %n",correct,errors,errorRate);

        return null;
    }

//    private TableReport makeTableReport(int[][] cnts) {
//        TableReport
//    }


    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Evaluate Imputation Accuracy";
    }

    @Override
    public String getToolTipText() {
        return "Evaluate Imputation Accuracy";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ImputationAccuracyPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }


}
