package me.ingiumdev.gobi.go;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class Mapping {
    Map<String, Set<Integer>> map;

    public Mapping() {
        this.map = new HashMap<>();
    }

    public Mapping(Map<String, Set<Integer>> map) {
        this.map = map;
    }

    public static Mapping createEnsemblMapping(String path) {
        Map<String, Set<Integer>> map = new HashMap<>();
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
                    Set<Integer> goTerms = map.get(hgnc);
                    if (goTerms == null) {
                        goTerms = new HashSet<>();
                        map.put(hgnc, goTerms);
                    }
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
                            int goTerm = Integer.parseInt(line.substring(i, k));
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
        Map<String, Set<Integer>> map = new HashMap<>();
        try (FileInputStream fileInputStream = new FileInputStream(path);
             GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
             BufferedReader br = new BufferedReader(new InputStreamReader(gzipInputStream))) {

            String line;
            String currentGene = null;

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

                // If we encounter a new gene, retrieve or create its list of GO terms
                if (currentGene == null || !currentGene.equals(geneName)) {
                    currentGene = geneName;
                    map.putIfAbsent(currentGene, new HashSet<>()); // Ensure a list exists for the gene
                }

                // Retrieve the list for the current gene
                Set<Integer> currentGOTerms = map.get(currentGene);

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
                        currentGOTerms.add(Integer.parseInt(goTermColumn.substring(goStart, goEnd)));
                        goStart = goEnd + 1; // Move past the '|' delimiter
                    }
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return new Mapping(map);
    }


    Set<Integer> getTerms(String hgnc) {
        return map.get(hgnc);
    }

    public Map<String, Set<Integer>> getMap() {
        return map;
    }
}
