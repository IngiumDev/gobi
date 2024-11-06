package gtf.structs;

import java.util.HashMap;
import java.util.Map;

public class GTFAttributes {
    public static final String GENE_ID = "gene_id";
    public static final String GENE_NAME = "gene_name";
    public static final String TRANSCRIPT_ID = "transcript_id";
    public static final String TRANSCRIPT_NAME = "transcript_name";
    public static final String EXON_ID = "exon_id";
    public static final String EXON_NUMBER = "exon_number";
    public static final String PROTEIN_ID = "protein_id";
    public static final String CCDS_ID = "ccds_id";
    public static final String CCDS_ID2 = "ccdsid";
    // Gene
    private String geneID;
    private String geneName;
    // Transcript
    private String transcriptID;
    private String transcriptName;
    // Exon
    private String exonID; //not always present
    private String exonNumber; //not always present
    // CDS
    private String proteinID; //not always present: rely on
    private String ccdsID; //not always present

    public GTFAttributes(Builder builder) {
        this.geneID = builder.geneID;
        this.geneName = builder.geneName;
        this.transcriptID = builder.transcriptID;
        this.transcriptName = builder.transcriptName;
        this.exonID = builder.exonID;
        this.exonNumber = builder.exonNumber;
        this.proteinID = builder.proteinID;
        this.ccdsID = builder.ccdsID;
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

    public static GTFAttributes parseAttributes(String attributeString) {
        // TODO: seperate by entry, if we have all required attributes, break;
        GTFAttributes.Builder attribute = new GTFAttributes.Builder();
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

            // Determine if the value is quoted
            String value;
            if (i < len && attributeString.charAt(i) == '"') {
                // The value should start with a quote
                int valueStart = ++i; // Skip the opening quote
                while (i < len && attributeString.charAt(i) != '"') {
                    i++;
                }
                if (i >= len) {
                    break; // Malformed input
                }
                value = attributeString.substring(valueStart, i++);
            } else {
                // Read until the next space or semicolon for unquoted values
                int valueStart = i;
                while (i < len && attributeString.charAt(i) != ' ' && attributeString.charAt(i) != ';') {
                    i++;
                }
                value = attributeString.substring(valueStart, i);
            }

            // Add the key-value pair to the map
            switch (key) {
                case GENE_ID -> attribute.setGeneID(value);
                case GENE_NAME -> attribute.setGeneName(value);
                case TRANSCRIPT_ID -> attribute.setTranscriptID(value);
                case TRANSCRIPT_NAME -> attribute.setTranscriptName(value);
                case EXON_ID -> attribute.setExonID(value);
                case EXON_NUMBER -> attribute.setExonNumber(value);
                case PROTEIN_ID -> attribute.setProteinID(value);
                case CCDS_ID, CCDS_ID2 -> attribute.setCcdsID(value);
            }

            // Skip any space after the value
            while (i < len && attributeString.charAt(i) == ' ') {
                i++;
            }

            // Skip the semicolon separating attributes
            if (i < len && attributeString.charAt(i) == ';') {
                i++;
            }
        }

        return attribute.build();
    }


    public String getGeneID() {
        return geneID;
    }

    public void setGeneID(String geneID) {
        this.geneID = geneID;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public void setTranscriptID(String transcriptID) {
        this.transcriptID = transcriptID;
    }

    public String getTranscriptName() {
        return transcriptName;
    }

    public void setTranscriptName(String transcriptName) {
        this.transcriptName = transcriptName;
    }

    public String getExonID() {
        return exonID;
    }

    public void setExonID(String exonID) {
        this.exonID = exonID;
    }

    public String getExonNumber() {
        return exonNumber;
    }

    public void setExonNumber(String exonNumber) {
        this.exonNumber = exonNumber;
    }

    public String getProteinID() {
        return proteinID;
    }

    public void setProteinID(String proteinID) {
        this.proteinID = proteinID;
    }

    public String getCcdsID() {
        return ccdsID;
    }

    public void setCcdsID(String ccdsID) {
        this.ccdsID = ccdsID;
    }

    static class Builder {
        // Gene
        private String geneID;
        private String geneName;
        // Transcript
        private String transcriptID;
        private String transcriptName;
        // Exon
        private String exonID; //not always present
        private String exonNumber; //not always present
        // CDS
        private String proteinID; //not always present: rely on
        private String ccdsID; //not always present

        public Builder() {
        }

        public void setGeneID(String geneID) {
            this.geneID = geneID;
        }

        public void setGeneName(String geneName) {
            this.geneName = geneName;
        }

        public void setTranscriptID(String transcriptID) {
            this.transcriptID = transcriptID;
        }

        public void setTranscriptName(String transcriptName) {
            this.transcriptName = transcriptName;
        }

        public void setExonID(String exonID) {
            this.exonID = exonID;
        }

        public void setExonNumber(String exonNumber) {
            this.exonNumber = exonNumber;
        }

        public void setProteinID(String proteinID) {
            this.proteinID = proteinID;
        }

        public void setCcdsID(String ccdsID) {
            this.ccdsID = ccdsID;
        }

        public GTFAttributes build() {
            return new GTFAttributes(this);
        }
    }
}
