package BaseComponents;

import ExonSkipping.ES_SE;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Genome {
    private final HashSet<Gene> allGenes;
    private final ArrayList<ES_SE> exonSpliceEvents;
    private final HashMap<String, Gene> id2Gene;
    private final HashMap<String, Transcript> id2transcript;

    public Genome() {
        allGenes = new HashSet<>();
        this.exonSpliceEvents = new ArrayList<>();
        this.id2Gene = new HashMap<>();
        this.id2transcript = new HashMap<>();
    }

    @Override
    public String toString(){
        return "Genome with " + this.getAllGenes().size() + "many genomes.";
    }

    // Getter and Setter
    public void addGene(Gene toAdd){
        this.allGenes.add(toAdd);
    }

    public HashSet<Gene> getAllGenes() {
        return allGenes;
    }

    public ArrayList<ES_SE> getExonSpliceEvents() {
        return exonSpliceEvents;
    }

    public HashMap<String, Gene> getId2Gene() {
        return id2Gene;
    }

    public HashMap<String, Transcript> getId2transcript() {
        return id2transcript;
    }
}
