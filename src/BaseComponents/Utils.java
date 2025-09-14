package BaseComponents;

import java.util.*;

/**
 * this class contains all methods, which can be used in different scenarios, exercises...
 */
public abstract class Utils {
    public static final String RESET =" \u001B[0m";
    public static final String RED = "\u001B[31m";
    public static final String GREEN ="\u001B[32m";
    public static final String WHITE ="\u001B[37m";
    public static final String YELLOW = "\u001B[33m";
    public static final String BLUE = "\u001B[34m";
    public static final String CYAN = "\u001B[36m";

    public final Map<Character, Character> complementMap = Map.of(
            'A', 'T',
            'T', 'A',
            'C', 'G',
            'G', 'C'
    );

    public static void checkIfNull(Object obj, Object hasToBe){
        if(obj == null){
            System.err.println("Element is null instead of " + hasToBe.getClass() + "!");
        }
    }

    public static String getReverseComplement(StringBuilder dna) {
        int length = dna.length();
        StringBuilder reverseComplement = new StringBuilder(length);

        // Traverse the DNA sequence in reverse order
        for (int i = length - 1; i >= 0; i--) {
            char base = dna.charAt(i);
            switch (base) {
                case '-':
                    reverseComplement.append('-');
                    break;
                case 'A':
                    reverseComplement.append('T');
                    break;
                case 'T':
                    reverseComplement.append('A');
                    break;
                case 'C':
                    reverseComplement.append('G');
                    break;
                case 'G':
                    reverseComplement.append('C');
                    break;
                default:
                    throw new IllegalArgumentException("Invalid DNA base: " + base);
            }
        }

        return reverseComplement.toString();
    }

    // brute force it since .equals()-method in Intervals does not work;
    public static boolean checkIfEqual(Collection<Interval> set1, Collection<Interval> set2) {
        A: for(Interval interval: set1) {
            for (Interval toComp : set2) {
                if (toComp.getStartGenomic() == interval.getStartGenomic() && toComp.getEndGenomic() == interval.getEndGenomic()) {
                    continue A;
                }
            }
            return false;
        }
        return true;
    }

    public static boolean checkIfContained(Set<TreeSet<Interval>> big, TreeSet<Interval> set2) {
        for(TreeSet<Interval> interval: big) {
            if(Utils.checkIfEqual(set2, interval)) return true;
        }
        return false;
    }

}
