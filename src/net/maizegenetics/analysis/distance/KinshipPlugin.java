package net.maizegenetics.analysis.distance;

import com.google.common.collect.Range;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
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
 * @author Terry Casstevens
 * @author Zhiwu Zhang
 * @author Peter Bradbury
 *
 * modified by Peter Bradbury, 9/26/2014 converted to self-describing Plugin and
 * to use Phenotype package
 */
public class KinshipPlugin extends AbstractPlugin {

    public static enum KINSHIP_METHOD {

        Scaled_IBS,
        GCTA
    };
    private GenotypeTable.GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GenotypeTable.GENOTYPE_TABLE_COMPONENT[]{
        GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT.ReferenceProbability, GenotypeTable.GENOTYPE_TABLE_COMPONENT.AlleleProbability};

    private PluginParameter<KINSHIP_METHOD> method = new PluginParameter.Builder<>("method", KINSHIP_METHOD.Scaled_IBS, KINSHIP_METHOD.class)
            .guiName("Kinship method")
            .range(Range.encloseAll(Arrays.asList(KINSHIP_METHOD.values())))
            .description("The Scaled_IBS (Endelman) method produces a kinship matrix that is scaled to give a reasonable estimate of additive genetic variance. "
                    + "The GCTA uses the algorithm published here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/pdf/main.pdf.")
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
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if ((alignInList == null) || (alignInList.isEmpty())) {
            throw new IllegalArgumentException("KinshipPlugin: Nothing selected. Please select a genotype.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

        List<Datum> result = new ArrayList<>();
        Iterator itr = alignInList.iterator();
        while (itr.hasNext()) {

            Datum current = (Datum) itr.next();
            String datasetName = current.getName();
            DistanceMatrix kin = null;

            if (current.getData() instanceof GenotypeTable) {
                GenotypeTable myGenotype = (GenotypeTable) current.getData();
                if (kinshipMethod() == KINSHIP_METHOD.Scaled_IBS) {
                    //kin = Kinship.createKinship(myGenotype, Kinship.KINSHIP_TYPE.Endelman, myDatatype.value());
                    kin = EndelmanDistanceMatrix.getInstance(myGenotype, 6, this);
                } else if (kinshipMethod() == KINSHIP_METHOD.GCTA) {
                    kin = GCTADistanceMatrix.getInstance(myGenotype, this);
                } else {
                    throw new IllegalArgumentException("Unknown method to calculate kinship: " + kinshipMethod());
                }
            } else {
                throw new IllegalArgumentException("Invalid selection. Can't create kinship matrix from: " + datasetName);
            }

            if (kin != null) {
                //add kin to datatree;
                Datum ds = new Datum("kin_" + datasetName, kin, "Kinship matrix created from " + datasetName);
                result.add(ds);
            }

        }

        return new DataSet(result, this);

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = KinshipPlugin.class.getResource("/net/maizegenetics/analysis/images/Kin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Kinship";
    }

    @Override
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
