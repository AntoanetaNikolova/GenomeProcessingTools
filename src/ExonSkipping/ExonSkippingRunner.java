package ExonSkipping;

import BaseComponents.Gene;
import BaseComponents.Intron;
import org.apache.commons.cli.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class ExonSkippingRunner {
    public static void main(String[] args) {
        Options options = new Options();
        options.addOption("gtf",true, "Input gtf file");
        options.addOption("o",  true, "Input output file path");
        CommandLineParser parser = new BasicParser();

        String gtfFilename = "";
        String outputFilename = "";

        try {
            CommandLine cmd = parser.parse(options, args);
            if(cmd.hasOption("gtf")) {
                gtfFilename = cmd.getOptionValue("gtf");
            }

            if(cmd.hasOption("o")){
                outputFilename = cmd.getOptionValue("o");
            }

            if(!(cmd.hasOption("o") || cmd.hasOption("gtf"))){
                System.out.println("-gtf <give GTF-file path>" + "\n" +
                        " -o <output-file path>");
                return;
            }
        } catch (ParseException e) {
            System.err.println("Error parsing command line arguments!");
            return;
        }

        long startTime = System.currentTimeMillis();
        ExonSkipping exonSkipping = new ExonSkipping(gtfFilename);
        String output = exonSkipping.defineAndOrganizeIntrons();

        long endTime = System.currentTimeMillis();
        long duration = endTime - startTime;
        System.out.println("Runtime: " + duration + " ms");

        try (FileWriter writer = new FileWriter(outputFilename)) {
            writer.write(output);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

}
