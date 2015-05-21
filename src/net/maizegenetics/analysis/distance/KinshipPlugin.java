package net.maizegenetics.analysis.distance;

import com.google.common.collect.Range;
import net.maizegenetics.analysis.distance.Kinship.KINSHIP_TYPE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.distance.DistanceMatrix;

import javax.swing.*;

import java.net.URL;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * Author: Zhiwu Zhang Date: Apr 29, 2007
 *
 * modified by Peter Bradbury, 9/26/2014 converted to self-describing Plugin and
 * to use Phenotype package
 */
public class KinshipPlugin extends AbstractPlugin {

    public static enum KINSHIP_METHOD {

        Scaled_IBS, Pairwise_IBS, Pedigree
    };
    private GenotypeTable.GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GenotypeTable.GENOTYPE_TABLE_COMPONENT[]{
        GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability};

    private PluginParameter<KINSHIP_METHOD> method = new PluginParameter.Builder<>("method", KINSHIP_METHOD.Scaled_IBS, KINSHIP_METHOD.class)
            .guiName("Kinship method")
            .range(Range.encloseAll(Arrays.asList(KINSHIP_METHOD.values())))
            .description("The scaled_IBS method produces a kinship matrix that is scaled to give a reasonable estimate of additive genetic variance. The pairwise_IBS method, which "
                    + "is the method used by TASSEL ver.4, may result in an inflated estimate of genetic variance. Either will do a good job of controlling population structure in MLM. "
                    + "The pedigree method is used to calculate a kinship matrix from a pedigree information.")
            .build();
    private PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> myDatatype = new PluginParameter.Builder<>("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.class)
            .genotypeTable()
            .range(GENOTYPE_COMP)
            .description("If the genotype table contains more than one type of genotype data, choose the type to use for calculating kinship.")
            .build();

    public KinshipPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input.getSize() == 0) {
            throw new IllegalArgumentException("KinshipPlugin: Nothing selected. Please select a genotype or pedigree data.");
        }
    }

    public DataSet processData(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataSet();

            List<Datum> result = new ArrayList<>();
            Iterator itr = alignInList.iterator();
            while (itr.hasNext()) {

                Datum current = (Datum) itr.next();
                String datasetName = current.getName();
                Kinship kin = null;

                try {

                    if (current.getData() instanceof GenotypeTable) {
                        //this section implements additional options for calculating kinship
                        GenotypeTable myGenotype = (GenotypeTable) current.getData();
                        if (method.value() == KINSHIP_METHOD.Pairwise_IBS) {
                            kin = new Kinship(myGenotype, KINSHIP_TYPE.IBS, myDatatype.value());
                        } else if (method.value() == KINSHIP_METHOD.Scaled_IBS) {
                            kin = new Kinship(myGenotype, KINSHIP_TYPE.Endelman, myDatatype.value());
                        } else {
                            throw new IllegalArgumentException("The pedigree method cannot be used to calculate kinship from genotype data.");
                        }
                    } else if (current.getData() instanceof Phenotype) { //pedigree data
                        Phenotype ped = (Phenotype) current.getData();
                        kin = new Kinship(ped);
                    } else {
                        String message = "Invalid selection. Can't create kinship matrix from: " + datasetName;
                        if (isInteractive()) {
                            JOptionPane.showMessageDialog(getParentFrame(), message);
                        } else {
                            System.out.println(message);
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    String message = "Problem creating kinship matrix from: " + datasetName + "\n" + e.getClass().getName() + ": " + e.getMessage();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), message);
                    } else {
                        System.out.println(message);
                        e.printStackTrace();
                    }
                }

                if (kin != null) {
                    //add kin to datatree;
                    Datum ds = new Datum("kin_" + datasetName, kin.getDm(), "kinship matrix created from " + datasetName);
                    result.add(ds);
                }

            }

            return new DataSet(result, this);

        } finally {
            fireProgress(100);
        }
    }

    public ImageIcon getIcon() {
        URL imageURL = KinshipPlugin.class.getResource("/net/maizegenetics/analysis/images/Kin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Kinship";
    }

    public String getToolTipText() {
        return "Calculate kinship from marker data";
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(KinshipPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    public DistanceMatrix runPlugin(DataSet input) {
        return (DistanceMatrix) performFunction(input).getData(0).getData();
    }

    /**
     * The scaled_IBS method produces a kinship matrix that is scaled to give a
     * reasonable estimate of additive genetic variance. The pairwise_IBS
     * method, which is the method used by TASSEL ver.4, may result in an
     * inflated estimate of genetic variance. Either will do a good job of
     * controlling population structure in MLM. The pedigree method is used to
     * calculate a kinship matrix from a pedigree information.
     *
     * @return Kinship method
     */
    public KINSHIP_METHOD kinshipMethod() {
        return method.value();
    }

    /**
     * Set Kinship method. The scaled_IBS method produces a kinship matrix that
     * is scaled to give a reasonable estimate of additive genetic variance. The
     * pairwise_IBS method, which is the method used by TASSEL ver.4, may result
     * in an inflated estimate of genetic variance. Either will do a good job of
     * controlling population structure in MLM. The pedigree method is used to
     * calculate a kinship matrix from a pedigree information.
     *
     * @param value Kinship method
     *
     * @return this plugin
     */
    public KinshipPlugin kinshipMethod(KINSHIP_METHOD value) {
        method = new PluginParameter<>(method, value);
        return this;
    }

    /**
     * If the genotype table contains more than one type of genotype data,
     * choose the type to use for calculating kinship.
     *
     * @return Genotype Component
     */
    public GENOTYPE_TABLE_COMPONENT genotypeComponent() {
        return myDatatype.value();
    }

    /**
     * Set Genotype Component. If the genotype table contains more than one type
     * of genotype data, choose the type to use for calculating kinship.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    public KinshipPlugin genotypeComponent(GENOTYPE_TABLE_COMPONENT value) {
        myDatatype = new PluginParameter<>(myDatatype, value);
        return this;
    }

}
