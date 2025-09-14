package BaseComponents;

import java.util.LinkedList;
import java.util.Random;

public class Read {
    // add all needed output information
    private final int readLength;

    // use for seq extraction and genomic region vectors
    private int fwStartRelative;
    private int fwEndRelative;
    private int rwStartRelative;
    private int rwEndRelative;

    // use for output
    private int fwTrStartRelative;
    private int fwTrEndRelative;
    private int rwTrStartRelative;
    private int rwTrEndRelative;

    private String fw;
    private String rw;
    private final LinkedList<Integer> mutationsfw;
    private final LinkedList<Integer> mutationsrw;

    private final Random random;
    private static final char[] bases = {'A', 'T', 'C', 'G'};


    public Read(int readLength){
        this.readLength = readLength;
        this.mutationsfw = new LinkedList<>();
        this.mutationsrw = new LinkedList<>();
        this.random = new Random();
        this.fwStartRelative = -1;
    }

    /**
     * @param dna input sequence to mutate
     * @param mutationRate at which rate to introduce a miss match
     * @param x on which strand the mutations are introduced: 'f' := fw and 'r' := rw
     * @return the sequence, but with introduced miss matches
     */
    public String mutateDNA(String dna, double mutationRate, char x) {
        StringBuilder mutatedDNA = new StringBuilder(dna.length());

        for (int i = 0; i < dna.length(); i++) {
            char base = dna.charAt(i);

            // Determine if the base should mutate
            if (random.nextInt(100) < mutationRate) {
                if(x == 'f') mutationsfw.add(i);
                else mutationsrw.add(i);
                // Mutate the base
                base = getRandomMutation(base, random);
            }

            mutatedDNA.append(base);
        }

        return mutatedDNA.toString();
    }

    /**
     * @param base original nucleotide
     * @param random object
     * @return new nucleotide after "mutation"
     */
    private static char getRandomMutation(char base, Random random) {
        // Generate a random base that's different from the original
        char newBase;
        do {
            newBase = bases[random.nextInt(bases.length)];
        } while (newBase == base);

        return newBase;
    }

    @Override
    public String toString(){
        return " ";
    }

    // Getter and Setter
    public String getFw() {
        return fw;
    }

    public String getRw() {
        return rw;
    }

    public void setFw(String fw, double mutationRate) {
        // introduce mutations and reverse complement
        fw = mutateDNA(fw,mutationRate,'f');
        this.fw = fw;
    }

    public void setRw(String rw, double mutationRate) {
        rw = mutateDNA(rw,mutationRate,'r');
        this.rw = rw;
    }

    public int getReadLength() {
        return readLength;
    }

    public int getFwStartRelative() {
        return fwStartRelative;
    }

    public LinkedList<Integer> getMutationsfw() {
        return mutationsfw;
    }

    public LinkedList<Integer> getMutationsrw() {
        return mutationsrw;
    }

    public int getRwStartRelative() {
        return rwStartRelative;
    }

    /** since all needed info is available, set indexes */
    public void setFragmentStartRelative(int fragmentStartRelative, int fragmentLength, int transcriptLength,  String strandDirection) {
        this.fwTrStartRelative = fragmentStartRelative;
        this.fwTrEndRelative = this.fwTrStartRelative + readLength;
        this.rwTrStartRelative = fragmentStartRelative + fragmentLength - readLength;
        this.rwTrEndRelative = fragmentStartRelative + fragmentLength;

        if (strandDirection.equals("+")) {
            this.fwStartRelative = fwTrStartRelative;
            this.fwEndRelative = fwTrEndRelative;
            this.rwStartRelative = rwTrStartRelative;
            this.rwEndRelative = rwTrEndRelative;
        } else {  // on the negative strand
            this.fwStartRelative = transcriptLength - this.fwTrEndRelative;
            this.fwEndRelative = this.fwStartRelative + readLength;
            this.rwStartRelative = transcriptLength - this.rwTrEndRelative;
            this.rwEndRelative = this.rwStartRelative + readLength;
        }
    }

    public int getFwEndRelative() {
        return fwEndRelative;
    }

    public int getRwEndRelative() {
        return rwEndRelative;
    }

    public int getFwTrEndRelative() {
        return fwTrEndRelative;
    }

    public int getFwTrStartRelative() {
        return fwTrStartRelative;
    }

    public int getRwTrEndRelative() {
        return rwTrEndRelative;
    }

    public int getRwTrStartRelative() {
        return rwTrStartRelative;
    }
}
