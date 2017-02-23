/*
 *  TaxaListUtils
 */
package net.maizegenetics.taxa;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;

/**
 *
 * @author Terry Casstevens
 */
public class TaxaListUtils {

    private TaxaListUtils() {
        // utility class
    }

    /**
     * Intersect joins the specified taxa.
     *
     * @param group1 an TaxaList
     * @param group2 another TaxaList
     *
     * @return the taxa in the intersection of groups 1 and 2, sorted in
     * ascending order
     */
    public static TaxaList getCommonTaxa(TaxaList group1, TaxaList group2) {
        return getCommonTaxa(new TaxaList[]{group1, group2});
    }

    /**
     * Intersect joins the specified taxa.
     *
     * @param groups groups to join.
     *
     * @return The taxa from the intersect join sorted alphabetically
     */
    public static TaxaList getCommonTaxa(TaxaList[] groups) {
        return getCommonTaxa(groups, true);
    }

    /**
     * Intersect joins the specified taxa.
     *
     * @param groups groups to join.
     * @param sorted whether to sort taxa alphabetically
     *
     * @return The taxa from the intersect join
     */
    public static TaxaList getCommonTaxa(TaxaList[] groups, boolean sorted) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        } else if (groups.length == 1 && !sorted) {
            return groups[0];
        }

        LinkedHashSet<Taxon> intersectIds = new LinkedHashSet<>();
        for (int x = 0; x < groups[0].numberOfTaxa(); x++) {
            intersectIds.add(groups[0].get(x));
        }
        for (int i = 1; i < groups.length; i++) {
            List<Taxon> temp = new ArrayList<>();
            for (int j = 0; j < groups[i].numberOfTaxa(); j++) {
                temp.add(groups[i].get(j));
            }
            intersectIds.retainAll(temp);
        }

        TaxaListBuilder builder = new TaxaListBuilder();
        builder.addAll(intersectIds);
        if (sorted) {
            builder.sortTaxaAlphabetically();
        }
        return builder.build();

    }

    /**
     * Union joins the specified taxa.
     *
     * @param group1 an TaxaList
     * @param group2 another TaxaList
     *
     * @return	the taxa in the union of taxa 1 and 2, sorted in ascending
     * order
     */
    public static TaxaList getAllTaxa(TaxaList group1, TaxaList group2) {
        return getAllTaxa(new TaxaList[]{group1, group2});
    }

    public static TaxaList getAllTaxa(TaxaList[] groups) {
        return getAllTaxa(groups, true);
    }

    /**
     * Union joins the specified taxa.
     *
     * @param groups taxa to join.
     * @param sorted whether to sort taxa alphabetically
     *
     * @return The taxa from the union join
     */
    public static TaxaList getAllTaxa(TaxaList[] groups, boolean sorted) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        }

        LinkedHashSet<Taxon> allIds = new LinkedHashSet<>();
        for (int i = 0; i < groups.length; i++) {
            int n = groups[i].numberOfTaxa();
            for (int j = 0; j < n; j++) {
                allIds.add(groups[i].get(j));
            }
        }

        TaxaListBuilder builder = new TaxaListBuilder();
        builder.addAll(allIds);
        if (sorted) {
            builder.sortTaxaAlphabetically();
        }
        return builder.build();

    }
}
