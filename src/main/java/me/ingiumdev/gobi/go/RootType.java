package me.ingiumdev.gobi.go;

public enum RootType {
    BIOLOGICAL_PROCESS('b'),
    CELLULAR_COMPONENT('c'),
    MOLECULAR_FUNCTION('m');

    private final char value;

    RootType(char value) {
        this.value = value;
    }

    public static RootType parseRootType(String text) {
        for (RootType type : RootType.values()) {
            if (type.value == text.charAt(0)) {
                return type;
            }
        }
        throw new IllegalArgumentException("Unknown RootType: " + text);
    }

    public char getValue() {
        return value;
    }
}
