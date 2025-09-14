package ReadSimulator;

import BaseComponents.Exon;
import java.util.Comparator;
import java.util.TreeSet;

public class PreTranscript {
    private final String chromosome;
    private final String strandDirection;
    private final String transcriptID;
    private final TreeSet<Exon> exons;
    private int transcriptLength;
    private boolean tooShort;

    public PreTranscript(String chr, String strandDirection, String geneID, String transcriptID){
        this.chromosome = chr;
        this.strandDirection = strandDirection;
        this.transcriptID = transcriptID;
        exons = new TreeSet<>(Comparator.comparingLong(Exon::getStartGenomic));
        this.tooShort = false;
    }


    // Getter and Setter
    public TreeSet<Exon> getExons() {
        return exons;
    }

    public void setTranscriptLength(int transcriptLength, int readLength) {
        this.transcriptLength = transcriptLength;
        if(transcriptLength < readLength) tooShort = true;
    }

    public long getTranscriptLength() {
        return transcriptLength;
    }

    public String getStrandDirection() {
        return strandDirection;
    }

    public String getChromosome() {
        return chromosome;
    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public boolean getTooShort(){
        return tooShort;
    }
}
