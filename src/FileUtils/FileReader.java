package FileUtils;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class FileReader {
    public static File getFile(String path) {
        return new File(path);
    }

    public static ArrayList<String> readLinesFromFile(String path) {
        return readLinesFromFile(getFile(path));
    }

    public static ArrayList<String> readLinesFromFile(File file) {
        ArrayList<String> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new java.io.FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                lines.add(line);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return lines;
    }

    public static ArrayList<String[]> readSeparatedLinesFromFile(File file, String separator) {
        ArrayList<String[]> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new java.io.FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                lines.add(line.split(separator));
            }
        } catch (Exception e) {
            e.printStackTrace();

        }
        return lines;
    }

    public static ArrayList<String[]> readCSVFromFile(File file) {
        return readSeparatedLinesFromFile(file, ",");
    }

    public static ArrayList<String[]> readTSVFromFile(File file) {
        return readSeparatedLinesFromFile(file, "\t");
    }

    public static Map<String, String> parseAttributes(String attributeString) {
        Map<String, String> attributes = new HashMap<>();
        String[] parts = attributeString.split(";");
        for (String part : parts) {
            part = part.trim();
            if (!part.isEmpty()) {
                String[] keyValue = part.split(" ", 2);
                if (keyValue.length == 2) {
                    String key = keyValue[0].trim();
                    String value = keyValue[1].trim().replaceAll("^\"|\"$", ""); // Remove surrounding quotes
                    attributes.put(key, value);
                }
            }
        }
        return attributes;
    }


}
