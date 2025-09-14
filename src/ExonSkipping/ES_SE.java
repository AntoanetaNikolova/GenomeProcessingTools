package ExonSkipping;

import BaseComponents.Intron;
import BaseComponents.Utils;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.stream.Collectors;

public class ES_SE {
    private String ID;  // gene_id
    private String symbol;  // gene name
    private String chromosome;  // chromosome
    private char strand;  // +/-
    private int n_prots;  // number of annotated CDS in the gene
    private int n_trans;  // number of annotated transcripts
    private ArrayList<Intron> SV;  // the SV-intron as "start:end"
    private ArrayList<Intron> WT;   // the WT introns within the SV-intron, separated by |
    private ArrayList<Intron> intronsInSV;
    private ArrayList<Intron> allIntronsInEvent;
    private ArrayList<String> SV_prots;  // ids of the SV CDS, separated by |
    private ArrayList<String> WT_prots;  // ids of the WT CDS, separated by |
    private int min_skipped_exon;
    private int max_skipped_exon;
    private int min_skipped_bases;
    private int max_skipped_bases;

    public ES_SE(String ID, String symbol, String chr, char strand, int n_prots, int n_trans, ArrayList<Intron> SV, ArrayList<Intron> WT, ArrayList<Intron> allIntronsInEvent){
        this.ID = ID;
        this.symbol = symbol;
        this.chromosome = chr;
        this.strand = strand;
        this.n_prots = n_prots;
        this.n_trans = n_trans;
        this.SV = SV;
        this.WT = WT;
        this.allIntronsInEvent = allIntronsInEvent;
        this.SV_prots = new ArrayList<>();
        this.WT_prots = new ArrayList<>();
        this.intronsInSV = new ArrayList<>();
    }


    @Override
    public String toString(){
        StringBuilder output = new StringBuilder(this.ID + "\t" + this.symbol + "\t" + this.chromosome + "\t" + this.strand + "\t" +
                this.n_prots + "\t" + this.n_trans + "\t" +
                (this.SV.getFirst().getStartGenomic() + 1) + ":" + this.SV.getFirst().getEndGenomic() + "\t" );

        ArrayList<Intron> uniqueIntervals = new ArrayList<>(
                intronsInSV.stream()
                        .collect(Collectors.toMap(
                                interval -> interval.getStartGenomic() + ":" + interval.getEndGenomic(),
                                interval -> interval,
                                (existing, replacement) -> existing
                        )).values());
        uniqueIntervals.sort(Comparator.comparingInt(Intron::getStartGenomic).reversed());

        for(Intron intron: allIntronsInEvent){
            output.append((intron.getStartGenomic()+1) + ":" + intron.getEndGenomic() + "|");
        }
        output.deleteCharAt(output.length()-1);
        output.append("\t");

        // SV_prots:
        for(Intron intron: WT){
            output.append(intron.getCds3End().getID() + "|");
        }
        output.deleteCharAt(output.length()-1);
        output.append("\t");

        for(Intron intron: SV){
            output.append(intron.getCds3End().getID() + "|");
        }
        output.deleteCharAt(output.length()-1);
        output.append("\t");

        // MIN/MAX Exon

        output.append(this.min_skipped_exon + "\t");
        output.append(this.max_skipped_exon + "\t");
        output.append((this.min_skipped_bases + 1) + "\t");
        output.append((this.max_skipped_bases + 1) + "\t");

        // MIN/MAX Bases


        return output.toString();
    }


    // Getter and Setter
    public int getMin_skipped_exon() {
        return min_skipped_exon;
    }

    public int getMax_skipped_exon() {
        return max_skipped_exon;
    }

    public int getMin_skipped_bases() {
        return min_skipped_bases;
    }

    public int getMax_skipped_bases() {
        return max_skipped_bases;
    }

    public ArrayList<Intron> getIntronsInSV() {
        return intronsInSV;
    }

    public void setIntronsInSV(ArrayList<Intron> intronsInSV) {
        this.intronsInSV = intronsInSV;
    }

    public void setMin_skipped_exon(int min_skipped_exon) {
        this.min_skipped_exon = min_skipped_exon;
    }

    public void setMax_skipped_exon(int max_skipped_exon) {
        this.max_skipped_exon = max_skipped_exon;
    }

    public void setMax_skipped_bases(int max_skipped_bases) {
        this.max_skipped_bases = max_skipped_bases;
    }

    public void setMin_skipped_bases(int min_skipped_bases) {
        this.min_skipped_bases = min_skipped_bases;
    }

    public String getID() {
        return ID;
    }
}
