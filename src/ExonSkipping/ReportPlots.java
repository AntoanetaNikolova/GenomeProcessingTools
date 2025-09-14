package ExonSkipping;

import BaseComponents.Intron;
import BaseComponents.Utils;

import java.io.IOException;
import java.util.*;
import java.io.FileWriter;
import java.util.stream.Collectors;

public class ReportPlots {
    public static void main(String[] args) {
        String[] GTFs = {
                "src/ExonSkipping/data/Homo_sapiens.GRCh37.67.gtf",
                "src/ExonSkipping/data/Homo_sapiens.GRCh37.75.gtf"
        };

        ArrayList<Integer>[] maxSkippedExons = new ArrayList[GTFs.length];
        for(int i = 0; i < GTFs.length; i++){
            maxSkippedExons[i] = new ArrayList<>();
        }
        ArrayList<Integer>[] maxSkippedBases = new ArrayList[GTFs.length];
        for(int i = 0; i < GTFs.length; i++){
            maxSkippedBases[i] = new ArrayList<>();
        }


        int gtfCounter = 0;
        ArrayList<ES_SE> allEvents = new ArrayList<>();

        for(String inputFilePath : GTFs) {
            // run all files and define ES-SE
            ExonSkipping exonSkipping = new ExonSkipping(inputFilePath);
            exonSkipping.defineAndOrganizeIntrons();

            for(ES_SE event: exonSkipping.getGenome().getExonSpliceEvents()){
                maxSkippedExons[gtfCounter].add(event.getMax_skipped_exon());
                maxSkippedBases[gtfCounter].add(event.getMax_skipped_bases());
            }

            gtfCounter++;

            // Genes test
            allEvents.addAll(exonSkipping.getGenome().getExonSpliceEvents());
        }


        // max skipped exons writer
        try (FileWriter writer = new FileWriter("src/ExonSkipping/data/exons.csv")) {
            for (ArrayList<Integer> valuesPerFile : maxSkippedExons) {
                StringBuilder line = new StringBuilder();
                for(Integer i: valuesPerFile){
                    line.append(i).append(",");
                }
                if(line.isEmpty()) writer.write("0");
                else line.deleteCharAt(line.length()-1);
                writer.write(line + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // max skipped bases writer
        try (FileWriter writer = new FileWriter("src/ExonSkipping/data/bases.csv")) {
            for (ArrayList<Integer> valuesPerFile : maxSkippedBases) {
                StringBuilder line = new StringBuilder();
                for(Integer i: valuesPerFile){
                    line.append(i).append(",");
                }
                if(line.isEmpty()) writer.write("0");
                else line.deleteCharAt(line.length()-1);
                writer.write(line + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // call python program on the saved values - line 1: plot for distribution of max skipped exons, line 2: -||- bases
        try {
            ProcessBuilder pb = new ProcessBuilder("python", "src/ExonSkipping/plot_max_exon_bases_data.py");
            pb.inheritIO();  // This will allow you to see any output or errors from the Python script
            Process process = pb.start();
            int exitCode = process.waitFor();  // Wait for the script to finish
            if (exitCode == 0) {
                System.out.println("Plot generated successfully.");
            } else {
                System.out.println("Error: Python script did not run successfully.");
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }

        // Determine top 10 genes
        // Step 1: map each GeneID to an Integer[]
        HashMap<String, Integer[]> geneId2numMaxExons = new HashMap<>();

        for(ES_SE event: allEvents){
            if(!geneId2numMaxExons.containsKey(event.getID())){
                geneId2numMaxExons.put(event.getID(), new Integer[]{event.getMax_skipped_exon(), event.getMin_skipped_bases()});
            } else {
                Integer[] oldValues = geneId2numMaxExons.get(event.getID());
                geneId2numMaxExons.get(event.getID())[0] = oldValues[0] + event.getMax_skipped_exon();
                geneId2numMaxExons.get(event.getID())[1] = oldValues[1] + event.getMin_skipped_bases();
            }
        }

        List<Map.Entry<String, Integer[]>> entryList = new ArrayList<>(geneId2numMaxExons.entrySet());

        // Step 2: Sort the list by the two integers in Integer[] in descending order
        entryList.sort((e1, e2) -> {
            Integer[] values1 = e1.getValue();
            Integer[] values2 = e2.getValue();

            // Compare by the first value in descending order
            int compareFirst = values2[0].compareTo(values1[0]);
            if (compareFirst != 0) {
                return compareFirst;
            }
            // If the first values are equal, compare by the second value in descending order
            return values2[1].compareTo(values1[1]);
        });

        // Step 3: Print the first 10 entries
        System.out.println("Top 10 genes:");
        entryList.stream()
                .limit(10)
                .forEach(entry -> System.out.println(entry.getKey() + " => " + Arrays.toString(entry.getValue())));

    }

}
