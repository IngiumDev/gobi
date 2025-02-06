package me.ingiumdev.gobi.go;

import java.util.*;

public class GOTerm {
    private final int id;
    private final String fullID;
    private final String name;
    private boolean isSOT; // SOT = Standard of Truth
    private final Set<String> associatedGenes;
    private final List<GOTerm> parents;
    private final List<GOTerm> children;
    private List<Integer> tempParentIDs;

    public String getFullID() {
        return fullID;
    }

    private GOTerm(Builder builder) {
        name = builder.name;
        id = builder.id;
        fullID = builder.fullID;
        tempParentIDs = builder.tempParentIDs;
        isSOT = false;
        parents = new ArrayList<>();
        children = new ArrayList<>();
        associatedGenes = new HashSet<>();
    }

    public void setSOT() {
        isSOT = true;
    }

    public boolean addGene(String gene) {
        return associatedGenes.add(gene);
    }

    public List<Integer> getTempParentIDs() {
        return tempParentIDs;
    }

    public List<GOTerm> getParents() {
        return parents;
    }

    public int getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public void addParent(GOTerm parent) {
        parents.add(parent);
    }

    public void addChild(GOTerm child) {
        children.add(child);
    }

    public boolean isSOT() {
        return isSOT;
    }

    public Set<String> getAssociatedGenes() {
        return associatedGenes;
    }

    public List<GOTerm> getChildren() {
        return children;
    }

    @Override
    public String toString() {
        return "GOEntry{" +
                "fullID='" + fullID + '\'' +
                ", name='" + name + '\'' +
                '}';
    }

    public void clearTempParentIDs() {
        tempParentIDs = null;
    }

    public static final class Builder {
        private String name;
        private String fullID;
        private int id;
        private List<Integer> tempParentIDs;

        public Builder() {
        }


        public Builder(GOTerm copy) {
            this.name = copy.getName();
            this.id = copy.getId();
            this.tempParentIDs = copy.getTempParentIDs();
        }

        public Builder setName(String name) {
            this.name = name;
            return this;
        }

        public Builder setId(String id) {
            this.id = Integer.parseInt(id);
            this.fullID = id;
            return this;
        }

        public Builder setTempParentIDs(List<Integer> tempParentIDs) {
            this.tempParentIDs = tempParentIDs;
            return this;
        }

        public GOTerm build() {
            return new GOTerm(this);
        }

        @Override
        public boolean equals(Object o) {
            if (o == null || getClass() != o.getClass()) return false;

            Builder builder = (Builder) o;
            return id == builder.id && Objects.equals(name, builder.name) && Objects.equals(fullID, builder.fullID) && Objects.equals(tempParentIDs, builder.tempParentIDs);
        }

        @Override
        public int hashCode() {
            return id;
        }
    }
}
