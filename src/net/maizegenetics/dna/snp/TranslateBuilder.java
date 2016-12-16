/*
 *  TranslateBuilder
 * 
 *  Created on Dec 10, 2016
 */
package net.maizegenetics.dna.snp;

/**
 *
 * @author Terry Casstevens
 */
public class TranslateBuilder {

    private TranslateBuilder() {
    }

    /**
     *
     * @param translateTaxa translate taxa
     * @param translateSite translate site
     *
     * @return Translate
     */
    public static Translate getInstance(TranslateIndex translateTaxa, TranslateIndex translateSite) {

        if (translateTaxa == null) {
            throw new IllegalArgumentException("TranslateBuilder: filter: must specify translateTaxa");
        }

        if (translateSite == null) {
            throw new IllegalArgumentException("TranslateBuilder: filter: must specify translateSite");
        }

        return new Translate(translateTaxa, translateSite);

    }

}
