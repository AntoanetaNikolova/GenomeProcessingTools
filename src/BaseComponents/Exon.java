package BaseComponents;

public class Exon {
    private final int startGenomic;
    private final int endGenomic;
    private int startRelative;
    private int endRelative;
    private final int length;

    public Exon(int startGenomic, int endGenomic) {
        this.startGenomic = startGenomic;
        this.endGenomic = endGenomic;
        this.length = endGenomic - startGenomic;
    }

    @Override
    public String toString() {
        return startRelative + ":" + endRelative + ", " + startGenomic + ":" + endGenomic;
    }

    // Getter and Setter
    public int getStartGenomic() {
        return startGenomic;
    }

    public int getEndGenomic() {
        return endGenomic;
    }

    public int getEndRelative() {
        return endRelative;
    }

    public int getStartRelative() {
        return startRelative;
    }

    public void setStartRelative(int startRelative) {
        this.startRelative = startRelative;
    }

    public void setEndRelative(int endRelative) {
        this.endRelative = endRelative;
    }

    public int getLength() {
        return length;
    }
}
