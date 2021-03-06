/*
 *  FilterByTaxa
 * 
 *  Created on Dec 21, 2016
 */
package net.maizegenetics.analysis.filter;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.FilterTaxa;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.ListStats;
import net.maizegenetics.dna.snp.genotypecall.Stats;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

/**
 *
 * @author Terry Casstevens
 */
public class FilterByTaxa {

    private FilterByTaxa() {
    }

    public static GenotypeTable filter(GenotypeTable orig, FilterTaxa filter) {

        int numSites = orig.numberOfSites();
        int numTaxa = orig.numberOfTaxa();
        TaxaList taxa = orig.taxa();

        List<String> taxaNamesFromFilter = filter.taxaNames();
        HashSet<String> taxaNames = null;
        if (taxaNamesFromFilter != null) {
            taxaNames = new HashSet<>(taxaNamesFromFilter);
        }

        TaxaList taxaList = filter.taxaList();
        if (taxaList != null) {
            if (taxaNames == null) {
                taxaNames = new HashSet<>();
            }
            for (Taxon current : taxaList) {
                taxaNames.add(current.getName());
            }
        }

        List<Integer> includedTaxa = null;
        if (taxaNames != null) {
            includedTaxa = new ArrayList<>();
            if (filter.includeTaxa()) {
                for (int t = 0; t < numTaxa; t++) {
                    Taxon taxon = taxa.get(t);
                    if (taxaNames.remove(taxon.getName())) {
                        includedTaxa.add(t);
                    }
                }
            } else {
                for (int t = 0; t < numTaxa; t++) {
                    Taxon taxon = taxa.get(t);
                    if (!taxaNames.remove(taxon.getName())) {
                        includedTaxa.add(t);
                    }
                }
            }
        }

        Stream<Stats> stream = null;

        if (filter.minNotMissing() != 0.0) {
            stream = stream(includedTaxa, orig, stream);
            stream = stream.filter((Stats stats) -> {
                return stats.percentNotMissing() >= filter.minNotMissing();
            });
        }

        if (filter.minHeterozygous() != 0.0 || filter.maxHeterozygous() != 1.0) {
            stream = stream(includedTaxa, orig, stream);
            stream = stream.filter((Stats stats) -> {
                double hetFreq = stats.proportionHeterozygous();
                return filter.minHeterozygous() <= hetFreq && filter.maxHeterozygous() >= hetFreq;
            });
        }

        TaxaList keepTaxa = null;
        if (stream == null) {
            TaxaListBuilder builder = new TaxaListBuilder();
            if (includedTaxa != null) {
                for (int current : includedTaxa) {
                    builder.add(taxa.get(current));
                }
            }
            keepTaxa = builder.build();
        } else {
            keepTaxa = stream.map((Stats stats) -> taxa.get(stats.index()))
                    .collect(TaxaList.collect());
        }

        if (keepTaxa.numberOfTaxa() == 0) {
            return null;
        } else if (keepTaxa.numberOfTaxa() == taxa.numberOfTaxa()) {
            return orig;
        } else {
            return FilterGenotypeTable.getInstance(orig, keepTaxa, false);
        }

    }

    private static Stream<Stats> stream(List<Integer> includedTaxa, GenotypeTable orig, Stream<Stats> stream) {
        if (stream != null) {
            return stream;
        }
        ListStats taxaStats = ListStats.getTaxaInstance(orig.genotypeMatrix());
        if (includedTaxa == null) {
            return IntStream.range(0, orig.numberOfTaxa()).parallel().mapToObj((int value) -> taxaStats.get(value));
        } else {
            return includedTaxa.stream().parallel().map((Integer t) -> taxaStats.get(t));
        }
    }
}
