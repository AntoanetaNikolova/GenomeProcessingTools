package BaseComponents;

public class Intron {
    private Transcript originTranscript;
    private final CDS cds5End;
    private final CDS cds3End;
    private final int length;
    private final int startGenomic;
    private final int endGenomic;
    private boolean checked = false;


    public Intron(CDS cds5end, CDS cds3End){
        this.cds3End = cds3End;
        this.cds5End = cds5end;
        this.endGenomic = cds3End.getStartGenomic();
        this.startGenomic = cds5end.getEndGenomic();
        this.length = endGenomic - startGenomic + 1;
    }

    @Override
    public String toString(){
        return "Intron on " + this.cds5End.getEndGenomic() + ":" + this.cds3End.getStartGenomic()
                + " from transcript " + this.originTranscript.getID();
    }


    // Setter and Getter
    public void setOriginTranscript(Transcript originTranscript) {
        this.originTranscript = originTranscript;
    }

    public Transcript getOriginTranscript() {
        return originTranscript;
    }

    public CDS getCds3End() {
        return cds3End;
    }

    public int getLength() {
        return length;
    }

    public void setChecked() {
        this.checked = true;
    }

    public boolean getChecked() {
        return this.checked;
    }

    public int getEndGenomic() {
        return endGenomic;
    }

    public int getStartGenomic() {
        return startGenomic;
    }
}
