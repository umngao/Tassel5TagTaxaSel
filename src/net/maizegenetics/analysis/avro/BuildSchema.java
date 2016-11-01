/*
 *  BuildSchema
 * 
 *  Created on Aug 31, 2016
 */
package net.maizegenetics.analysis.avro;

import java.io.File;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.LoggingUtils;
import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;

/**
 *
 * @author Terry Casstevens
 */
public class BuildSchema {

    public static void main(String[] args) {

        LoggingUtils.setupDebugLogging();

        try {

            GenotypeTable genotype = ImportUtils.read("mdp_genotype.hmp.txt");
            int numTaxa = genotype.numberOfTaxa();
            int numSites = genotype.numberOfSites();

            SchemaBuilder.FieldAssembler<Schema> genotypeSchemaBuilder = SchemaBuilder
                    .builder("net.maizegenetics")
                    .record("genotype")
                    .fields();
            for (int s = 0; s < numSites; s += AvroConstants.GENOTYPE_BLOCK_SIZE) {
                for (int t = 0; t < numTaxa; t += AvroConstants.GENOTYPE_BLOCK_SIZE) {
                    genotypeSchemaBuilder = genotypeSchemaBuilder
                            .name(AvroConstants.getKey(t, s))
                            .type(AvroConstants.BYTE_BLOCK_SCHEMA)
                            .noDefault();
                }
            }
            Schema genotypeSchema = genotypeSchemaBuilder.endRecord();

            Schema tasselSchema = SchemaBuilder
                    .builder("net.maizegenetics")
                    .record("tassel")
                    .fields()
                    .name("taxa").type(AvroConstants.TAXA_SCHEMA).noDefault()
                    .name("positions").type(AvroConstants.POSITIONS_SCHEMA).noDefault()
                    .name("genotype").type(genotypeSchema).noDefault()
                    .endRecord();

            System.out.println("Genotype Schema...");
            System.out.println(genotypeSchema.toString(true));
            System.out.println("Tassel Schema...");
            System.out.println(tasselSchema.toString(true));
            System.out.println("schema name: " + tasselSchema.getName());
            System.out.println("schema type: " + tasselSchema.getType());
            //System.out.println("schema namespace: " + tasselSchema.getNamespace());
            System.out.println("schema fullname: " + tasselSchema.getFullName());

            DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(tasselSchema);
            try (DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<>(datumWriter)) {
                dataFileWriter.setCodec(CodecFactory.snappyCodec());
                dataFileWriter.create(tasselSchema, new File("DS_1.avro"));
                dataFileWriter.append(new GenericRecordGenotypeTable(tasselSchema, genotype));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
}
