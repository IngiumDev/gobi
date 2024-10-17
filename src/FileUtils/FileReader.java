package FileUtils;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;

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
}
