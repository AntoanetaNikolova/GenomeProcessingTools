package BaseComponents;

import augmentedTree.DefaultInterval;

public class Interval implements augmentedTree.Interval {
    private final int startGenomic;
    private int endGenomic;
    private final int length;

    public Interval(int start, int end) {
        this.startGenomic = start;
        this.endGenomic = end;
        this.length = startGenomic - endGenomic;
    }

    @Override
    public String toString(){
        return startGenomic + "-" + endGenomic;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof DefaultInterval))
            return false;
        @SuppressWarnings("unchecked")
        Interval o = (Interval)obj;
        return o.startGenomic==startGenomic && o.endGenomic==endGenomic;
    }

    public int hashCode() {
        return startGenomic+(endGenomic<<13);
    }

    public int getLength() {
        return length;
    }

    public int getEndGenomic() {
        return endGenomic;
    }

    public int getStartGenomic() {
        return startGenomic;
    }

    @Override
    public int getStart() {
        return startGenomic;
    }

    @Override
    public int getStop() {
        return endGenomic;
    }

    public void setStop(int stop) {
        this.endGenomic = stop;
    }
}
