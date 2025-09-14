package BaseComponents;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

public class Transcript {
    private Gene originGene;
    private String ID;
    private String name;
    private TreeSet<CDS> allCDS;
    private ArrayList<Intron> allIntrons;
    private ArrayList<Integer> intronStartingPositions;
    private final TreeSet<Exon> exons;
    private long startGenomic;
    private long startRelative;
    private long endGenomic;
    private long endRelative;
    private long chromosomeStart;

    public Transcript() {
        this.allCDS = new TreeSet<>(Comparator.comparingLong(CDS::getStartGenomic));
        this.allIntrons = new ArrayList<>();
        this.exons = new TreeSet<>(Comparator.comparingLong(Exon::getStartGenomic));
        this.startRelative = Integer.MAX_VALUE;
        this.endRelative = Integer.MIN_VALUE;
    }

    public Transcript(String id, String name){
        this();
        this.ID = id;
        this.name = name;
    }

    public Transcript(String ID){
        this.ID = ID;
        this.exons = new TreeSet<>(Comparator.comparingLong(Exon::getStartGenomic).thenComparing(Exon::getEndGenomic));
    }

    @Override
    public String toString(){
        return "Transcript " + this.getID() + ", " + this.getName() + " with " + this.allCDS.size() + " many CDS.";
    }


    // Getter and Setter
    public String getID() {
        return ID;
    }

    public String getName() {
        return name;
    }

    public void addCDS(CDS toAdd){
        this.allCDS.add(toAdd);
    }

    public TreeSet<CDS> getAllCDS() {
        return allCDS;
    }

    public ArrayList<Intron> getAllIntrons() {
        return allIntrons;
    }

    public ArrayList<Integer> getIntronStartingPositions() {
        return intronStartingPositions;
    }

    public void setIntronStartingPositions(ArrayList<Integer> intronStartingPositions) {
        this.intronStartingPositions = intronStartingPositions;
    }

    public TreeSet<Exon> getExons() {
        return exons;
    }

    public long getStartGenomic() {
        return startGenomic;
    }

    public long getEndGenomic() {
        return endGenomic;
    }

    public void setIndexesGenomic(){
        this.startGenomic = exons.getFirst().getStartGenomic();
        this.endGenomic = exons.getLast().getEndGenomic();
    }

    public TreeSet<Interval> getAlignedExons(boolean strandPositive, int readStart, int readEnd){
        ArrayList<Exon> orderedExons = new ArrayList<>(exons);
        TreeSet<Interval> subIntervals = new TreeSet<>(Comparator.comparingInt(Interval::getStart).thenComparingInt(Interval::getStop));

        Exon exon;
        for(int i = 0; i < orderedExons.size(); i++){
            if(strandPositive) {
                exon = orderedExons.get(i);
            } else exon = orderedExons.get(orderedExons.size() - i - 1);

            if (readStart >= exon.getStartGenomic() && readStart <= exon.getEndGenomic() && readEnd >= exon.getStartGenomic() && readEnd <= exon.getEndGenomic()) {
                subIntervals.add(new Interval(readStart, readEnd));
                break;
            }
            else if (readStart >= exon.getStartGenomic() && readStart <= exon.getEndGenomic()) {
                subIntervals.add(new Interval(readStart, exon.getEndGenomic()));
            }
            else if (readStart <= exon.getStartGenomic() && readEnd >= exon.getEndGenomic()) {
                subIntervals.add(new Interval(exon.getStartGenomic(), exon.getEndGenomic()));
            }
            else if (readEnd >= exon.getStartGenomic() && readEnd <= exon.getEndGenomic()) {
                subIntervals.add(new Interval(exon.getStartGenomic(), readEnd));
                break;
            }
        }

        return subIntervals;
    }
}
