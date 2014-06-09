package net.maizegenetics.dna.tag;

import com.google.common.collect.Multiset;

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 *
 * @author Ed Buckler
 */
public interface TaxaDistribution {

    TaxaDistribution increment(int taxaNum);

    int[] depths();

    int[][] taxaWithDepths();

    int[] encodeTaxaDepth();

    Multiset<Integer> taxaDepthMap();

    int totalDepth();

    int maxTaxa();
}
