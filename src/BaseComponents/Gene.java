package BaseComponents;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Gene {
    private String ID;
    private String Name;
    private String chromosome;
    private String strandDirection;
    private final ArrayList<Transcript> allTranscripts;
    private ArrayList<Intron> allIntrons;
    private Map<Integer, ArrayList<Intron>> organizedIntrons;
    private int n_ports;
    private int n_trans;

    public Gene() {
        this.allTranscripts = new ArrayList<>();
        this.allIntrons = new ArrayList<>();
        this.organizedIntrons = new HashMap<>();
    }

    public Gene(String chromosome, String strand){
        this();
        this.chromosome = chromosome;
        this.strandDirection = strand;
    }


    @Override
    public String toString(){
        return "Gene with " + this.getID() + ", " + this.getName() + " with " +
                this .getAllTranscripts().size() + " many transcripts.";
    }

    // Getter and Setter
    public String getID() {
        return ID;
    }

    public void setID(String ID) {
        this.ID = ID;
    }

    public String getName() {
        return Name;
    }

    public void setName(String name) {
        Name = name;
    }

    public void addTranscript(Transcript toAdd){
        this.allTranscripts.add(toAdd);
    }

    public ArrayList<Transcript> getAllTranscripts() {
        return allTranscripts;
    }

    public ArrayList<Intron> getAllIntrons() {
        return allIntrons;
    }

    public void setAllIntrons(ArrayList<Intron> allIntrons) {
        this.allIntrons = allIntrons;
    }

    public void setOrganizedIntrons(Map<Integer, ArrayList<Intron>> organizedIntrons) {
        this.organizedIntrons = organizedIntrons;
    }

    public Map<Integer, ArrayList<Intron>> getOrganizedIntrons() {
        return organizedIntrons;
    }

    public String getChromosome() {
        return chromosome;
    }

    public String getStrandDirection() {
        return strandDirection;
    }

    public int getN_ports() {
        return n_ports;
    }

    public void setN_ports(int n_ports) {
        this.n_ports = n_ports;
    }

    public void setN_trans(int n_trans) {
        this.n_trans = n_trans;
    }

    public int getN_trans() {
        return n_trans;
    }
}
