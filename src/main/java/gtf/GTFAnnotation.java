package gtf;

import gtf.structs.Gene;

import java.util.HashMap;
import java.util.Map;

public class GTFAnnotation {
    private Map<String, Gene> genes = new HashMap<>();

    public Map<String, Gene> getGenes() {
        return genes;
    }

    public void setGenes(HashMap<String, Gene> genes) {
        this.genes = genes;
    }

    public void addGene(Gene gene) {
        genes.put(gene.getGeneID(), gene);
    }

}
