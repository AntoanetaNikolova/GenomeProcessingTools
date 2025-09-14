package BaseComponents;

public class CDS {
    private String ID;
    private final int startGenomic;   // []
    private final int endGenomic;
    private final int length;

    public CDS( int startGenomic, int endGenomic){
        this.endGenomic = endGenomic;
        this.startGenomic = startGenomic;
        this.length = endGenomic - startGenomic;
    }

    @Override
    public String toString(){
        return "CDS " + this.ID + " with " + this.startGenomic + ":" + this.endGenomic;
    }


    // Getter and Setter
    public void setID(String ID) {
        this.ID = ID;
    }

    public String getID() {
        return ID;
    }

    public int getEndGenomic() {
        return endGenomic;
    }

    public int getStartGenomic() {
        return startGenomic;
    }

    public int getLength() {
        return length;
    }
}
