package me.ingiumdev.gobi.parsers;

import me.ingiumdev.gobi.go.GOTerm;
import me.ingiumdev.gobi.go.RootType;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class OboParser {
    BufferedReader br;
    RootType rootType;
    private final HashMap<Integer, GOTerm> entries;

    public OboParser(String path, RootType rootType) {
        this.rootType = rootType;
        try {
            br = new BufferedReader(new FileReader(path));
            skipToNextTerm();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        entries = new HashMap<>();
    }

    private boolean skipToNextTerm() throws IOException {
        String line;
        while ((line = br.readLine()) != null) {
            if (!line.isEmpty()) {
                if (line.charAt(0) == '[' && line.charAt(2) == 'e') {
                    break;
                }
            }
        }
        return line == null;
    }

    public HashMap<Integer, GOTerm> parse() {
        boolean hasNext = processEntry();
        while (hasNext) {
            hasNext = processEntry();
        }
        return entries;
    }

    private boolean processEntry() {
        try {
            String test = br.readLine();
            String rawID = test.substring(7);
            String name = br.readLine().substring(6);
            if (doesNameSpaceMatch()) {
                List<Integer> parentIDs = new ArrayList<>();
                String line;
                while ((line = br.readLine()) != null) {
                    if (!line.isEmpty()) {
                        if (line.charAt(0) == '[') {
                            break;
                        } else if (line.charAt(0) == 'i') {
                            switch (line.charAt(3)) {
                                case 'a' ->
                                        parentIDs.add(Integer.parseInt(line.split(" ")[1].substring(3))); // TODO: optimize this
                                case 'o' -> {
                                    return checkForRemainingTerms();
                                }
                            }
                        }
                    }
                }
                GOTerm entry = new GOTerm.Builder().setId(rawID).setName(name).setTempParentIDs(parentIDs).build();
                entries.put(Integer.parseInt(rawID), entry);
                if (line == null) {
                    br.close();
                    return false;
                } else {
                    return true;
                }
            } else {
                return checkForRemainingTerms();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private boolean checkForRemainingTerms() throws IOException {
        boolean isOver = skipToNextTerm();
        if (isOver) {
            br.close();
            return false;
        } else {
            return processEntry();
        }
    }

    private boolean doesNameSpaceMatch() throws IOException {
        return br.readLine().charAt(11) == rootType.getValue();
    }
}