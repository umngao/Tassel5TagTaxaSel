/*
 *  FilterTaxaBuilderPlugin
 * 
 *  Created on Nov 29, 2016
 */
package net.maizegenetics.analysis.filter;

import com.google.common.collect.Range;
import java.awt.Frame;
import java.lang.reflect.Field;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.FilterTaxa;
import net.maizegenetics.dna.snp.FilterTaxa.FILTER_TAXA_ATTRIBUTES;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.TaxaList;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class FilterTaxaBuilderPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterTaxaBuilderPlugin.class);

    private PluginParameter<String> myFilterName = new PluginParameter.Builder<>(FILTER_TAXA_ATTRIBUTES.filterName.name(), "Filter", String.class)
            .description("Filter Name")
            .build();
    private PluginParameter<Double> myMinNotMissing = new PluginParameter.Builder<Double>(FILTER_TAXA_ATTRIBUTES.minNotMissing.name(), TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT, Double.class)
            .guiName("Min Proportion of Sites Present").range(Range.closed(0.0, 1.0)).build();
    private PluginParameter<Double> myMinHeterozygous = new PluginParameter.Builder<Double>(FILTER_TAXA_ATTRIBUTES.minHeterozygous.name(), TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT, Double.class)
            .guiName("Min Heterozygous Proportion").range(Range.closed(0.0, 1.0)).build();
    private PluginParameter<Double> myMaxHeterozygous = new PluginParameter.Builder<Double>(FILTER_TAXA_ATTRIBUTES.maxHeterozygous.name(), TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT, Double.class)
            .guiName("Max Heterozygous Proportion").range(Range.closed(0.0, 1.0)).build();
    private PluginParameter<Boolean> myIncludeTaxa = new PluginParameter.Builder<>(FILTER_TAXA_ATTRIBUTES.includeTaxa.name(), true, Boolean.class)
            .description("Include taxa from list of names or taxa list if true. Exclude otherwise.")
            .build();
    private PluginParameter<TaxaList> myTaxaList = new PluginParameter.Builder<>(FILTER_TAXA_ATTRIBUTES.taxaList.name(), null, TaxaList.class)
            .taxaList()
            .description("Filter based on taxa list.")
            .build();
    private PluginParameter<String> myTaxaNamesList = new PluginParameter.Builder<>(FILTER_TAXA_ATTRIBUTES.taxaNames.name(), null, String.class)
            .taxaNameList()
            .dependentOnParameter(myTaxaList, TAXA_LIST_NONE)
            .description("Filter based on taxa names.")
            .build();

    public FilterTaxaBuilderPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public FilterTaxaBuilderPlugin() {
        this(null, false);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        if (input == null) {
            return;
        }

    }

    @Override
    public DataSet processData(DataSet input) {

        Map<String, Object> values = new LinkedHashMap<>();
        for (Field field : getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(this);
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }
            if (current != null) {
                if (((current.value() != null)) && (!current.value().equals(current.defaultValue()))
                        || (current.cmdLineName().equals(FilterTaxa.FILTER_TAXA_ATTRIBUTES.filterName.name()))) {
                    values.put(current.cmdLineName(), current.value());
                }
            }
        }

        Object sites = values.get(FilterTaxa.FILTER_TAXA_ATTRIBUTES.taxaNames.name());
        if (sites != null) {
            values.put(FilterTaxa.FILTER_TAXA_ATTRIBUTES.taxaNames.name(), getListFromString(sites.toString()));
        }

        List<Datum> result = new ArrayList<>();

        FilterTaxa filter = new FilterTaxa(values);

        List<Datum> genotypeTableList = input.getDataOfType(GenotypeTable.class);
        if (genotypeTableList.size() >= 1) {
            for (Datum datum : genotypeTableList) {
                GenotypeTable current = (GenotypeTable) datum.getData();
                GenotypeTable filteredGenotype = GenotypeTableUtils.filter(filter, current);
                if ((filteredGenotype == null) || filteredGenotype.numberOfTaxa() == 0) {
                    DialogUtils.showWarning("No genotype data remained after filtering: " + datum.getName(), getParentFrame());
                } else if (filteredGenotype != current) {
                    Datum temp = new Datum(datum.getName() + "_" + filter.filterName(), filteredGenotype, null);
                    result.add(temp);
                    GenotypeSummaryPlugin.printSimpleSummary(temp);
                } else {
                    result.add(datum);
                    DialogUtils.showWarning("Genotype data unchanged after filtering: " + datum.getName(), getParentFrame());
                }
            }
        }

        List<Datum> phenoGenoTableList = input.getDataOfType(GenotypePhenotype.class);
        if (phenoGenoTableList.size() >= 1) {
            for (Datum datum : phenoGenoTableList) {
                GenotypePhenotype pheno = (GenotypePhenotype) datum.getData();
                GenotypeTable current = pheno.genotypeTable();
                GenotypeTable filteredGenotype = GenotypeTableUtils.filter(filter, current);
                if ((filteredGenotype == null) || filteredGenotype.numberOfTaxa() == 0) {
                    DialogUtils.showWarning("No genotype data remained after filtering: " + datum.getName(), getParentFrame());
                } else if (filteredGenotype != current) {
                    GenotypePhenotype resultPheno = new GenotypePhenotypeBuilder()
                            .genotype(filteredGenotype)
                            .phenotype(pheno.phenotype())
                            .union()
                            .build();
                    String name = datum.getName() + "_" + filter.filterName();
                    Datum temp = new Datum(name, resultPheno, null);
                    result.add(temp);
                    GenotypeSummaryPlugin.printSimpleSummary(filteredGenotype, name);
                } else {
                    result.add(datum);
                    DialogUtils.showWarning("Genotype data unchanged after filtering: " + datum.getName(), getParentFrame());
                }
            }
        }

        result.add(new Datum(filter.filterName(), new FilterList(filter), null));

        return new DataSet(result, this);

    }

    private List<String> getListFromString(String str) {

        if ((str == null) || (str.length() == 0)) {
            return null;
        }
        String[] tokens = str.split(",");
        List<String> result = new ArrayList<>();
        for (String current : tokens) {
            current = current.trim();
            if (current.length() != 0) {
                result.add(current);
            }
        }
        return result;

    }

    private String getStringFromList(List<String> list) {

        if ((list == null) || (list.isEmpty())) {
            return null;
        }

        StringBuilder builder = new StringBuilder();
        boolean first = true;
        for (String current : list) {
            current = current.trim();
            if (current.length() != 0) {
                if (first) {
                    first = false;
                } else {
                    builder.append(",");
                }
                builder.append(current);
            }
        }
        return builder.toString();

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = FilterTaxaBuilderPlugin.class.getResource("/net/maizegenetics/analysis/images/FilterNew.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Filter Genotype Table Taxa";
    }

    @Override
    public String getToolTipText() {
        return "Filter Genotype Table Taxa";
    }

    public FilterTaxa build() {
        return (FilterTaxa) performFunction(null).getData(0).getData();
    }

    public FilterTaxaBuilderPlugin useFilterValues(FilterTaxa filter) {
        setParametersToDefault();
        filter.attributes().entrySet().stream().forEach((attribute) -> {
            setParameter(attribute.getKey().name(), attribute.getValue());
        });
        return this;
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public FilterTaxa runPlugin() {
        return (FilterTaxa) performFunction(null).getData(0).getData();
    }

    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getDataOfType(GenotypeTable.class).get(0).getData();
    }

    public GenotypeTable runPlugin(GenotypeTable input) {
        return (GenotypeTable) performFunction(DataSet.getDataSet(input)).getDataOfType(GenotypeTable.class).get(0).getData();
    }

    /**
     * Filter Name
     *
     * @return Filter Name
     */
    public String filterName() {
        return myFilterName.value();
    }

    /**
     * Set Filter Name. Filter Name
     *
     * @param value Filter Name
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin filterName(String value) {
        myFilterName = new PluginParameter<>(myFilterName, value);
        return this;
    }

    /**
     * Min Proportion of Sites Present
     *
     * @return Min Proportion of Sites Present
     */
    public Double minNotMissing() {
        return myMinNotMissing.value();
    }

    /**
     * Set Min Proportion of Sites Present. Min Proportion of Sites Present
     *
     * @param value Min Proportion of Sites Present
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin minNotMissing(Double value) {
        myMinNotMissing = new PluginParameter<>(myMinNotMissing, value);
        return this;
    }

    /**
     * Min Heterozygous Proportion
     *
     * @return Min Heterozygous Proportion
     */
    public Double minHeterozygous() {
        return myMinHeterozygous.value();
    }

    /**
     * Set Min Heterozygous Proportion. Min Heterozygous Proportion
     *
     * @param value Min Heterozygous Proportion
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin minHeterozygous(Double value) {
        myMinHeterozygous = new PluginParameter<>(myMinHeterozygous, value);
        return this;
    }

    /**
     * Max Heterozygous Proportion
     *
     * @return Max Heterozygous Proportion
     */
    public Double maxHeterozygous() {
        return myMaxHeterozygous.value();
    }

    /**
     * Set Max Heterozygous Proportion. Max Heterozygous Proportion
     *
     * @param value Max Heterozygous Proportion
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin maxHeterozygous(Double value) {
        myMaxHeterozygous = new PluginParameter<>(myMaxHeterozygous, value);
        return this;
    }

    /**
     * Include Taxa
     *
     * @return Include Taxa
     */
    public Boolean includeTaxa() {
        return myIncludeTaxa.value();
    }

    /**
     * Set Include Taxa.
     *
     * @param value Include Taxa
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin includeTaxa(Boolean value) {
        myIncludeTaxa = new PluginParameter<>(myIncludeTaxa, value);
        return this;
    }

    /**
     * Taxa List
     *
     * @return Position List
     */
    public TaxaList taxaList() {
        return myTaxaList.value();
    }

    /**
     * Set Taxa List.
     *
     * @param value Taxa List
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin taxaList(TaxaList value) {
        myTaxaList = new PluginParameter<>(myTaxaList, value);
        return this;
    }

    /**
     * Taxa Names List
     *
     * @return Taxa Names List
     */
    public String taxaNamesList() {
        return myTaxaNamesList.value();
    }

    /**
     * Set Taxa Names List. Taxa Names List
     *
     * @param value Taxa Names List
     *
     * @return this plugin
     */
    public FilterTaxaBuilderPlugin taxaNamesList(String value) {
        myTaxaNamesList = new PluginParameter<>(myTaxaNamesList, value);
        return this;
    }
}
