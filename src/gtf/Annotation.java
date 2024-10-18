package gtf;

import java.util.HashMap;

public class Annotation {
    private HashMap<String, Gene> genes = new HashMap<String, Gene>();
    public HashMap<String, Gene> getGenes() {
        return genes;
    }

    public void setGenes(HashMap<String, Gene> genes) {
        this.genes = genes;
    }
    public void addGene(Gene gene) {
        genes.put(gene.getId(), gene);
    }

}
