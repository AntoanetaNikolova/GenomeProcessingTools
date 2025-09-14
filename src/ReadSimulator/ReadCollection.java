package ReadSimulator;

public class ReadCollection {
    private String chromosome;
    private String originGeneID;
    private String originTranscriptID;
    private int numReads;

    public ReadCollection(){}

    public ReadCollection(String gene, String originTranscript, int reads) {
        this.originGeneID = gene;
        this.originTranscriptID = originTranscript;
        this.numReads = reads;
    }

    // Getter and Setter
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getChromosome() {
        return chromosome;
    }

    public String getOriginGeneID() {
        return originGeneID;
    }

    public String getOriginTranscriptID() {
        return originTranscriptID;
    }

    public int getNumReads() {
        return numReads;
    }
}