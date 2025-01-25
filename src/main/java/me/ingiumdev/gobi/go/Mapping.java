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

    public Mapping CreateEnsemblMapping(String path) {
        Mapping mapping = new Mapping();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                int i = 0;
                int length = line.length();
                while (i < length) {
                    // iterate through the first gene ID
                    while (i < length && line.charAt(i) != '\t') {
                        i++;
                    }
                    i++; // skip to first char of hgnc
                    //  // Find the end of the hgnc
                    j = i;
                    while (i < length && line.charAt(i) != '\t') {
                        i++;
                    }
                    if (j == i) {
                        break;  // Empty column
                    }
                    String hgnc = line.substring(j, i++); // in first char of gos
                    List<String> terms = new ArrayList<>();


                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e)
        }

    }

    List<String> getTerms(String hgnc) {
        return map.get(hgnc);
    }
}
