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
    @Deprecated
    public static Map<String, String> parseAttributesWithRegex(String attributeString) {
        Map<String, String> attributes = new HashMap<>();
        if (attributeString.endsWith(";")) {
            attributeString = attributeString.substring(0, attributeString.length() - 1);
        }
        String[] parts = attributeString.trim().split("; ");

        for (String part : parts) {
            if (!part.isEmpty()) {
                String[] keyValue = part.split(" ");
                if (keyValue.length == 2) {
                    String key = keyValue[0];
                    String value = keyValue[1].substring(1, keyValue[1].length() - 1); // Remove surrounding quotes
                    attributes.put(key, value);
                }
            }
        }
        return attributes;
    }

    public static Map<String, String> parseAttributes(String attributeString) {
        Map<String, String> attributes = new HashMap<>();
        int len = attributeString.length();

        // Trim trailing semicolon, if it exists
        if (len > 0 && attributeString.charAt(len - 1) == ';') {
            attributeString = attributeString.substring(0, len - 1);
            len--;
        }
        int i = 0;
        while (i < len) {
            // Skip any leading spaces
            while (i < len && attributeString.charAt(i) == ' ') {
                i++;
            }

            // Find the end of the key
            int start = i;
            while (i < len && attributeString.charAt(i) != ' ') {
                i++;
            }
            if (start == i) {
                break; // No more key-value pairs
            }
            String key = attributeString.substring(start, i);

            // Skip space between key and value
            while (i < len && attributeString.charAt(i) == ' ') {
                i++;
            }

            // The value should start with a quote
            if (i >= len || attributeString.charAt(i) != '"') {
                break; // Malformed input
            }

            // Find the end of the quoted value
            int valueStart = ++i; // Skip the opening quote
            while (i < len && attributeString.charAt(i) != '"') {
                i++;
            }
            if (i >= len) {
                break; // Malformed input
            }
            String value = attributeString.substring(valueStart, i++);

            // Add the key-value pair to the map
            attributes.put(key, value);

            // Skip any space after the closing quote
            while (i < len && attributeString.charAt(i) == ' ') {
                i++;
            }

            // Skip the semicolon separating attributes
            if (i < len && attributeString.charAt(i) == ';') {
                i++;
            }
        }

        return attributes;
    }


}
