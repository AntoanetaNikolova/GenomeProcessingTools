package ReadSimulator;

import org.apache.commons.cli.*;

public class ReadSimulatorPlotsRunner {
    // read input
    public static void main(String[] args) throws Exception {
        Options options = new Options();
        options.addOption("length",true, "Input read length");
        options.addOption("frlength", true, "Input mean of normal distribution for fragment length");
        options.addOption("SD", true, "Input standard deviation of normal distribution for fragment length");
        options.addOption("mutationrate", true, "Input rate, at which bases in the sequence are mutated");
        options.addOption("readcounts", true, "Input file path for number of reads per transcript");
        options.addOption("fasta", true, "Input FASTA (.fa) file path");
        options.addOption("fidx", true, "Input FASTA index (.fai) file path");
        options.addOption("gtf", true, "Input GTF (.gtf) file path");
        options.addOption("od", true, "Input output file path");
        CommandLineParser parser = new BasicParser();

        int readLength = 0;
        double meanReadLength = 0.0;
        double SDReadLength = 0.0;
        String readCountsFilePath = "";
        double mutationsrate = 0.0;
        String fastaFilePath = "";
        String fidxFilePath = "";
        String gtfFilePath = "";
        String outputFilePath = "";

        try {
            CommandLine cmd = parser.parse(options, args);

            if(cmd.hasOption("length")) {
                String lengthInput = cmd.getOptionValue("length");
                try {
                    readLength = Integer.parseInt(lengthInput);
                } catch (NumberFormatException e) {
                    System.err.println("The option -length requires an integer. Invalid input: " + lengthInput);
                    System.exit(1);
                }
            }

            if(cmd.hasOption("frlength")) {
                String frlengthInput = cmd.getOptionValue("frlength");
                try {
                    meanReadLength = Double.parseDouble(frlengthInput);
                } catch (NumberFormatException e) {
                    System.err.println("The option -length requires an double. Invalid input: " + frlengthInput);
                    System.exit(1);
                }
            }

            if(cmd.hasOption("SD")) {
                String sdInput = cmd.getOptionValue("SD");
                try {
                    SDReadLength = Double.parseDouble(sdInput);
                } catch (NumberFormatException e) {
                    System.err.println("The option -length requires a double. Invalid input: " + sdInput);
                    System.exit(1);
                }
            }

            if(cmd.hasOption("readcounts")) {
                readCountsFilePath = cmd.getOptionValue("readcounts");
            }

            if(cmd.hasOption("mutationrate")) {
                String mutationsrateString = cmd.getOptionValue("mutationrate");
                try {
                    mutationsrate = Double.parseDouble(mutationsrateString);
                } catch (NumberFormatException e) {
                    System.err.println("The option -length requires a double. Invalid input: " + mutationsrateString);
                    System.exit(1);
                }
            }

            if(cmd.hasOption("fasta")) {
                fastaFilePath = cmd.getOptionValue("fasta");
            }

            if(cmd.hasOption("fidx")) {
                fidxFilePath = cmd.getOptionValue("fidx");
            }

            if(cmd.hasOption("gtf")){
                gtfFilePath = cmd.getOptionValue("gtf");
            }

            if(cmd.hasOption("od")){
                outputFilePath = cmd.getOptionValue("od");
            }

            if(!(cmd.hasOption("length") || cmd.hasOption("frlength") || cmd.hasOption("SD")
                    || cmd.hasOption("readcounts") || cmd.hasOption("mutationsrate") || cmd.hasOption("fasta")
                    || cmd.hasOption("fidx") || cmd.hasOption("gtf") || cmd.hasOption("od"))){

                System.out.println("Please, check all input arguments!\n" +
                        "-length <integer, read length>\n" +
                        "-frlength <double, mean normal distribution>\n" +
                        "-SD <double, standard deviation of normal distribution>\n" +
                        "-mutationrate <double, rate of mutated bases>}\n" +
                        "-readcounts <input file path for number of reads per transcript>\n" +
                        "-fasta <fasta input file path>\n" +
                        "-fidx <fasta input file path>\n" +
                        "-gtf <give GTF-file path>\n" +
                        "-o <output-file path>");
                System.exit(1);
            }
        } catch (ParseException e) {
            System.err.println("Error parsing command line arguments!");
            return;
        }


        ReadSimulatorPlots readSimulator = new ReadSimulatorPlots(readLength, meanReadLength, SDReadLength, readCountsFilePath,
                mutationsrate, fastaFilePath, fidxFilePath, gtfFilePath, outputFilePath);


        readSimulator.defineReadCollections();
        readSimulator.definePreTranscriptsFromGTF();
        readSimulator.modifyPreTranscript();
        readSimulator.createPlotsFrLength();
        readSimulator.createPlotsMutations();
        readSimulator.createPlotsReadStatistics();
        readSimulator.createPlotsMutationPosition();
    }
}

