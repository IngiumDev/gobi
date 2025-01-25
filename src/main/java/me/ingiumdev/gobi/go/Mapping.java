package me.ingiumdev.gobi.go;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class Mapping {
    Map<String, List<String>> map;

    public Mapping() {
        this.map = new HashMap<>();
    }

    public Mapping(Map<String, List<String>> map) {
        this.map = map;
    }

    public static Mapping createEnsemblMapping(String path) {
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

    public static Mapping createGOMapping(String path) {
        Map<String, List<String>> map = new HashMap<>();
        try (FileInputStream fileInputStream = new FileInputStream(path);
             GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
             BufferedReader br = new BufferedReader(new InputStreamReader(gzipInputStream))) {

            String line;
            String currentGene = null;
            List<String> currentGOTerms = null;

            while ((line = br.readLine()) != null) {
                int i = 0;
                int length = line.length();
                String geneName = null, column4 = null, goTermColumn = null;

                // Parse columns manually using iteration
                for (int col = 1; i < length; col++) {
                    int start = i;

                    // Move to the next tab or end of line
                    while (i < length && line.charAt(i) != '\t') {
                        i++;
                    }
                    String value = line.substring(start, i);

                    // Assign relevant columns based on column index
                    if (col == 3) {
                        geneName = value; // Column 3: Gene Name
                    } else if (col == 4) {
                        column4 = value; // Column 4: Skip condition
                    } else if (col == 5) {
                        goTermColumn = value; // Column 5: GO term
                    }

                    i++; // Skip the tab
                }

                // Skip line if column 4 is not empty
                if (column4 != null && !column4.isEmpty()) {
                    continue;
                }

                // If we encounter a new gene, save the previous one
                if (currentGene == null || !currentGene.equals(geneName)) {
                    if (currentGene != null && currentGOTerms != null) {
                        map.put(currentGene, currentGOTerms);
                    }
                    currentGene = geneName;
                    currentGOTerms = new ArrayList<>();
                }

                // Extract GO term numbers from the GO column
                if (goTermColumn != null && goTermColumn.startsWith("GO:")) {
                    int goStart = 3; // Skip "GO:"
                    while (goStart < goTermColumn.length()) {
                        int goEnd = goStart;

                        // Find the end of the GO term (either '|' or end of string)
                        while (goEnd < goTermColumn.length() && goTermColumn.charAt(goEnd) != '|') {
                            goEnd++;
                        }

                        // Extract and add the GO number
                        currentGOTerms.add(goTermColumn.substring(goStart, goEnd));
                        goStart = goEnd + 1; // Move past the '|' delimiter
                    }
                }
            }

            // Add the last gene processed
            if (currentGene != null && currentGOTerms != null) {
                map.put(currentGene, currentGOTerms);
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return new Mapping(map);
    }

    public static void main(String[] args) {
//        Mapping map = createEnsemblMapping("~/IdeaProjects/gobi/data/GOEnrich/goa_human_ensembl.tsv");
        Mapping map = createGOMapping("~/IdeaProjects/gobi/data/GOEnrich/goa_human.gaf.gz");
        System.out.println();
    }

    List<String> getTerms(String hgnc) {
        return map.get(hgnc);
    }
}
