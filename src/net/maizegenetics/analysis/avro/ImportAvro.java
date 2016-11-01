/*
 *  ImportAvro
 * 
 *  Created on Oct 11, 2016
 */
package net.maizegenetics.analysis.avro;

import java.util.Iterator;
import java.util.Map;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.genotypecall.GOBIIAvroGenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;
import org.apache.avro.generic.GenericArray;
import org.apache.avro.generic.GenericRecord;

/**
 *
 * @author Terry Casstevens
 */
public class ImportAvro {

    private ImportAvro() {
        // utility
    }

    public static GenotypeTable genotypeTable(GenericRecord genotypeTable) {

        GenotypeTable temp = ImportUtils.read("mdp_genotype.hmp.txt");
        PositionList positions = temp.positions();
        TaxaList taxa = taxa((GenericArray<GenericRecord>) genotypeTable.get(AvroConstants.GENOTYPE_TABLE_COMPONENTS.taxa.name()));
        GenotypeCallTable genotypes = genotypeCallTable(taxa.numberOfTaxa(), positions.numberOfSites(), false, (GenericRecord) genotypeTable.get(AvroConstants.GENOTYPE_TABLE_COMPONENTS.genotype.name()));
        return GenotypeTableBuilder.getInstance(genotypes, positions, taxa, null, null, null, null, null);

    }

    private static GenotypeCallTable genotypeCallTable(int numTaxa, int numSites, boolean phased, GenericRecord genotypes) {
        return GOBIIAvroGenotypeCallTable.getInstance(numTaxa, numSites, phased, genotypes);
    }

    private static TaxaList taxa(GenericArray<GenericRecord> taxa) {

        TaxaListBuilder builder = new TaxaListBuilder();
        Iterator<GenericRecord> itr = taxa.iterator();
        while (itr.hasNext()) {
            GenericRecord current = itr.next();
            builder.add(new Taxon(current.get(AvroConstants.TAXON_INDICES.name.name()).toString(), annotations((Map<String, String>) current.get(AvroConstants.TAXON_INDICES.annotations.name()))));
        }
        return builder.build();

    }

    private static GeneralAnnotation annotations(Map<String, String> annotations) {

        if (annotations.isEmpty()) {
            return null;
        }

        GeneralAnnotationStorage.Builder builder = GeneralAnnotationStorage.getBuilder();
        for (Map.Entry<String, String> current : annotations.entrySet()) {
            builder.addAnnotation(current.getKey(), current.getValue());
        }
        return builder.build();

    }

}
