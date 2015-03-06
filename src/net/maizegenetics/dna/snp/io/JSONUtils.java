/*
 *  JSONUtils
 * 
 *  Created on Mar 6, 2015
 */
package net.maizegenetics.dna.snp.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import javax.json.Json;
import javax.json.JsonArray;
import javax.json.JsonObject;
import javax.json.JsonReader;
import javax.json.JsonString;
import javax.json.stream.JsonGenerator;
import javax.json.stream.JsonGeneratorFactory;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class JSONUtils {

    private static final Logger myLogger = Logger.getLogger(JSONUtils.class);

    private JSONUtils() {
        // utility
    }

    private static void taxaListToJSON(TaxaList taxa, JsonGenerator generator) {
        generator.writeStartArray("TaxaList");
        taxa.stream().forEach((current) -> {
            taxaToJSON(current, generator);
        });
        generator.writeEnd();
    }

    private static void taxaToJSON(Taxon taxon, JsonGenerator generator) {
        generator.writeStartObject();
        generator.write("name", taxon.getName());
        generalAnnotationToJSON(taxon.getAnnotation(), generator);
        generator.writeEnd();
    }

    private static void generalAnnotationToJSON(GeneralAnnotation annotation, JsonGenerator generator) {
        Set<String> keys = annotation.getAnnotationKeys();
        if (keys.isEmpty()) {
            return;
        }
        generator.writeStartObject("anno");
        keys.stream().forEach((key) -> {
            String[] values = annotation.getTextAnnotation(key);
            if (values.length == 1) {
                generator.write(key, values[0]);
            } else {
                generator.writeStartArray(key);
                for (String current : values) {
                    generator.write(current);
                }
                generator.writeEnd();
            }
        });
        generator.writeEnd();
    }

    /**
     * Exports given taxa list to JSON file.
     *
     * @param taxa taxa list
     * @param filename file name (adds .json if needed)
     *
     * @return final filename
     */
    public static String exportTaxaListToJSON(TaxaList taxa, String filename) {
        filename = Utils.addSuffixIfNeeded(filename, ".json");
        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            Map<String, Object> properties = new HashMap<>(1);
            properties.put(JsonGenerator.PRETTY_PRINTING, true);
            JsonGeneratorFactory factory = Json.createGeneratorFactory(properties);
            try (JsonGenerator generator = factory.createGenerator(writer)) {
                generator.writeStartObject();
                taxaListToJSON(taxa, generator);
                generator.writeEnd();
            }
            return filename;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TaxaListIOUtils: exportTaxaListToJSON: problem saving file: " + filename + "\n" + e.getMessage());
        }
    }

    /**
     * Imports taxa list from JSON file.
     *
     * @param filename filename
     *
     * @return taxa list
     */
    public static TaxaList importTaxaListFromJSON(String filename) {
        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            try (JsonReader jsonReader = Json.createReader(reader)) {
                JsonObject obj = jsonReader.readObject();
                JsonArray taxaList = obj.getJsonArray("TaxaList");
                if (taxaList == null) {
                    throw new IllegalArgumentException("TaxaListIOUtils: importTaxaListFromJSON: There is no TaxaList in this file: " + filename);
                }
                return taxaListFromJSON(taxaList);
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TaxaListIOUtils: importTaxaListFromJSON: problem reading file: " + filename + "\n" + e.getMessage());
        }
    }

    private static TaxaList taxaListFromJSON(JsonArray taxaList) {
        if (taxaList.isEmpty()) {
            return null;
        }
        TaxaListBuilder builder = new TaxaListBuilder();
        taxaList.forEach((current) -> {
            builder.add(taxonFromJSON((JsonObject) current));
        });
        return builder.build();
    }

    private static Taxon taxonFromJSON(JsonObject json) {
        String name = json.getString("name");
        if (name == null) {
            throw new IllegalStateException("TaxaListIOUtils: taxonFromJSON: All Taxa must have a name.");
        }

        JsonObject anno = json.getJsonObject("anno");
        if (anno != null) {
            return new Taxon(name, generalAnnotationFromJSON(anno));
        } else {
            return new Taxon(name);
        }
    }

    private static GeneralAnnotation generalAnnotationFromJSON(JsonObject json) {
        GeneralAnnotationStorage.Builder builder = GeneralAnnotationStorage.getBuilder();
        json.forEach((key, value) -> {
            if (value instanceof JsonArray) {
                JsonArray jsonArray = (JsonArray) value;
                jsonArray.forEach((str) -> {
                    builder.addAnnotation(key, ((JsonString) str).getString());
                });
            } else if (value instanceof JsonString) {
                builder.addAnnotation(key, ((JsonString) value).getString());
            } else {
                throw new IllegalArgumentException("TaxaListIOUtils: generalAnnotationFromJSON: unknown value type: " + value.getClass().getName());
            }
        });
        return builder.build();
    }

}
