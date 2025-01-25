package me.ingiumdev.gobi.go;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Mapping {
    Map<String, List<String>> map;

    public Mapping() {
        this.map = new HashMap<>();
    }

    public Mapping(Map<String, List<String>> map) {
        this.map = map;
    }

    public static Mapping CreateEnsemblMapping(String path) {
        Map<String, List<String>> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line = br.readLine(); // Read the header or first line if needed
            while ((line = br.readLine()) != null) {
                int i = 0;
                int length = line.length();

                // Parse the first gene ID and hgnc columns
                while (i < length) {
                    // Iterate to the first tab (gene ID column)
                    while (i < length && line.charAt(i) != '\t') {
                        i++;
                    }
                    i++; // Skip to the first character of hgnc

                    // Find the end of the hgnc column
                    int j = i;
                    while (i < length && line.charAt(i) != '\t') {
                        i++;
                    }

                    if (j == i) {
                        break; // Empty column
                    }

                    String hgnc = line.substring(j, i++); // Move to first char of GO terms column

                    // Parse the GO terms
                    List<String> goTerms = new ArrayList<>();
                    while (i < length) {
                        // Look for the GO prefix
                        if (line.startsWith("GO:", i)) {
                            i += 3; // Skip past "GO:"
                            int k = i;

                            // Extract the GO term number until delimiter or end of column
                            while (k < length && line.charAt(k) != '|' && line.charAt(k) != '\t') {
                                k++;
                            }

                            // Add the GO term to the list
                            String goTerm = line.substring(i, k);
                            goTerms.add(goTerm);

                            i = k; // Move index forward
                        } else {
                            i++; // Skip non-GO text
                        }

                        // If a delimiter is found, skip past it
                        if (i < length && line.charAt(i) == '|') {
                            i++;
                        }
                    }

                    // Print or process the parsed results
                    map.put(hgnc, goTerms);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return new Mapping(map);
    }

    List<String> getTerms(String hgnc) {
        return map.get(hgnc);
    }

    public static void main(String[] args) {
        Mapping map = CreateEnsemblMapping("~/IdeaProjects/gobi/data/GOEnrich/goa_human_ensembl.tsv");
        System.out.println();
    }
}
