package ReadSimulator;

import BaseComponents.Exon;
import BaseComponents.Interval;
import BaseComponents.Read;
import BaseComponents.Utils;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import java.io.File;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.distribution.NormalDistribution;

public class ReadSimulatorPlots {
    // command line inputs
    private final String readCountsFilePath;
    private final String fastaFilePath;
    private final String findxFilePath;
    private final String gtfFilePath;
    private final int readLength;
    private final double mutationRate;
    private final String outputPath;

    // definition of ReadCollection objects
    private final HashMap<String, ReadCollection> transcriptID2ReadCollection;
    private final HashMap<String, PreTranscript> transcriptId2preTranscript;

    private final NormalDistribution distribution;
    private final FASTAIndex FASTAIndexFile;

    // Save data for plots
    ArrayList<Long> fragmentLengths;
    ArrayList<Integer> numberOfMutations;
    ArrayList<Integer> mutationPosition;
    long numAllReads;
    long numNONSplit;
    long numNONSplitNOMm;
    long numSplit;
    long numSplitNOMm;
    long numSplitNOMsRegions;


    public ReadSimulatorPlots(int length, double frlength, double SD, String readcounts , double mutationRate, String fasta, String fixd, String gtf, String outputPath) {
        this.readLength = length;
        this.readCountsFilePath = readcounts;
        this.mutationRate = mutationRate;
        this.fastaFilePath = fasta;
        this.findxFilePath = fixd;
        this.gtfFilePath = gtf;
        this.outputPath = outputPath;

        this.transcriptID2ReadCollection = new HashMap<>();
        this.transcriptId2preTranscript = new HashMap<>();
        this.distribution = new NormalDistribution(frlength, SD);
        FASTAIndexFile = new FASTAIndex(findxFilePath);

        // plot attributes
        this.fragmentLengths = new ArrayList<>();
        this.numberOfMutations = new ArrayList<>();
        this.numAllReads = 0;
        this.mutationPosition = new ArrayList<>();
    }

    /** create a ReadCollection object for each line in readcounts-file */
    protected void defineReadCollections() {
        try (BufferedReader br = new BufferedReader(new FileReader(readCountsFilePath))) {
            String line;
            while ((line = br.readLine()) != null){
                String[] elements = line.split("\t");
                if(elements[2].equals("count")) continue;
                if (elements.length != 3) System.err.println("Error! Line in readcounts-file with less/more than 3 elements");

                ReadCollection collection = new ReadCollection(elements[0], elements[1], Integer.parseInt(elements[2]));
                transcriptID2ReadCollection.put(elements[1], collection);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /** extract all needed information from the GTF, needed for read extraction */
    protected void definePreTranscriptsFromGTF(){
        try (BufferedReader br = new BufferedReader(new FileReader(gtfFilePath))) {
            String line;
            while ((line = br.readLine()) != null){
                int firstTab = line.indexOf("\t");
                if (firstTab != -1) {
                    if(!FASTAIndexFile.getChromosomes().contains(line.substring(0, firstTab).trim())) continue;
                    int secondTab = line.indexOf("\t", firstTab + 1);
                    if (secondTab != -1) {
                        String afterSecondTab = line.substring(secondTab + 1);

                        if (afterSecondTab.startsWith("exon")){
                            String[] values = line.split("\t");
                            String[] attributes = values[8].split("; ");

                            String currentTranscriptID = null;
                            String currentGeneID = null;

                            for(String attr: attributes){
                                if(attr.startsWith("gene_id")){
                                    currentGeneID = attr.substring(attr.indexOf(" ") + 1).replace("\"", "").replace(";", "");
                                } else if(attr.startsWith("transcript_id")){
                                    currentTranscriptID = attr.substring(attr.indexOf(" ") + 1).replace("\"", "").replace(";", "");
                                }

                                if(currentGeneID != null && currentTranscriptID != null) break;
                            }

                            if(!transcriptID2ReadCollection.containsKey(currentTranscriptID)) continue;
                            if(currentGeneID == null){
                                currentGeneID = attributes[0].split(" ")[2].replace("\"", "").replace(";", "");
                            }

                            if(!transcriptId2preTranscript.containsKey(currentTranscriptID)){ // if the transcript has not been defined => do it
                                PreTranscript preTranscript = new PreTranscript(values[0], values[6], currentGeneID, currentTranscriptID);
                                transcriptId2preTranscript.put(currentTranscriptID, preTranscript);
                            }
                            this.transcriptId2preTranscript.get(currentTranscriptID).getExons().add(new Exon(Integer.parseInt(values[3]), Integer.parseInt(values[4]) + 1));  // 1-based, end-exclusive
                        }
                    }
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Iteration through all transcripts and for each one:
     * 1. Define relative coordinated to each exon
     * 2. extract pre-transcript sequence from FASTA file using FASTA Index file
     * 3. modify pre-transcript sequence - remove introns
     * 4. get the corresponding ReadCollection object; for each read-object define the fragment and read indexes
     * 5. output
     **/
    protected void modifyPreTranscript() throws Exception {
        // Declare variables to avoid excessive memory allocation
        ArrayList<Interval> fwGenomicRegions;
        ArrayList<Interval> rwGenomicRegions;

        for(PreTranscript preTr: transcriptId2preTranscript.values()){
            // 1. Define relative coordinates exons
            int transcriptLength = 0;
            int startRelativeIndex = 0;

            for(Exon sortedExon: preTr.getExons()){
                sortedExon.setStartRelative(startRelativeIndex);
                int endRelative = startRelativeIndex + sortedExon.getLength();
                sortedExon.setEndRelative(endRelative);
                startRelativeIndex = endRelative;
                transcriptLength += sortedExon.getLength();
            }

            // If the length of the transcript is smaller than the read length => not possible to extract reads
            preTr.setTranscriptLength(transcriptLength, readLength);
            if(preTr.getTooShort()) {
                transcriptID2ReadCollection.remove(preTr.getTranscriptID());  // there is no need to iterate it, no reads will be created from it
                continue;
            }

            // 2. Extract pre-transcript sequence from FASTA using FASTA Index
            long preTranscriptStartGenomic = 0;
            long preTranscriptEndGenomic = 0;

            // 1-based and end-inclusive indexes
            preTranscriptStartGenomic = preTr.getExons().getFirst().getStartGenomic();
            preTranscriptEndGenomic = preTr.getExons().getLast().getEndGenomic() - 1;
            long preTranscriptLength = preTranscriptEndGenomic - preTranscriptStartGenomic + 1;

            String preTranscriptSequence = "";
            try {
                File fastaFile = new File(fastaFilePath);
                File faiFile = new File(findxFilePath);
                FastaSequenceIndex FASTAIndexFile = new FastaSequenceIndex(faiFile);
                IndexedFastaSequenceFile FASTAFile = new IndexedFastaSequenceFile(fastaFile, FASTAIndexFile);

                ReferenceSequence refSeq = FASTAFile.getSubsequenceAt(preTr.getChromosome(), preTranscriptStartGenomic, preTranscriptEndGenomic);

                // Output the sequence
                preTranscriptSequence = refSeq.getBaseString();
                FASTAFile.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            // 3. Modify pre-transcript sequence - remove introns
            StringBuilder preTranscriptStrBuilder = new StringBuilder(preTranscriptSequence);

            Exon end5 = null;
            int startPos = 0;
            int endPos;
            for(Exon end3: preTr.getExons()){
                if(end5 == null){
                    end5 = end3;
                    startPos = end3.getLength();
                    continue;
                }

                endPos = startPos + Math.abs(end3.getStartGenomic() - end5.getEndGenomic());
                preTranscriptStrBuilder.replace(startPos, endPos, "-".repeat(endPos-startPos));
                startPos = endPos + end3.getLength();
                end5 = end3;
            }

            preTranscriptSequence = preTranscriptStrBuilder.toString().replaceAll("-", "");
            if(preTranscriptSequence.length() != preTr.getTranscriptLength()) throw new Exception("Error while removing introns!");

            // 4. Get the corresponding ReadCollection object; for each read-object define the fragment and read indexes
            ReadCollection readColCurTranscript = transcriptID2ReadCollection.get(preTr.getTranscriptID());
            Utils.checkIfNull(readColCurTranscript, new ReadCollection());
            if(readColCurTranscript == null) continue;

            readColCurTranscript.setChromosome(preTr.getChromosome());

            Read r;
            for(int i = 0; i < readColCurTranscript.getNumReads(); i++){
                r = new Read(readLength);

                // Determine fragment length
                long fragmentLength = Math.round(distribution.sample());
                while(fragmentLength < readLength || fragmentLength > preTr.getTranscriptLength()) {
                    fragmentLength = Math.round(distribution.sample());
                }
                this.fragmentLengths.add(fragmentLength);

                // Determine fragment start
                if(fragmentLength == preTr.getTranscriptLength()) {
                    r.setFragmentStartRelative(0, (int)fragmentLength, (int)preTr.getTranscriptLength(), preTr.getStrandDirection());
                } else {
                    r.setFragmentStartRelative(Math.max(0, new Random().nextInt((int)(preTr.getTranscriptLength() - fragmentLength))), (int)fragmentLength, (int)preTr.getTranscriptLength(), preTr.getStrandDirection());
                }

                // Get genomic region vectors
                fwGenomicRegions = regionVectors(r.getFwStartRelative(), r.getFwEndRelative(), preTr);
                rwGenomicRegions = regionVectors(r.getRwStartRelative(), r.getRwEndRelative(), preTr);

                // Set Fw and Rw
                String leftEnd = preTranscriptSequence.substring(r.getFwStartRelative(), r.getFwEndRelative());
                String rightEnd = Utils.getReverseComplement(new StringBuilder(preTranscriptSequence.substring(r.getRwStartRelative(), r.getRwEndRelative())));

                if(preTr.getStrandDirection().equals("+")){
                    r.setFw(leftEnd, mutationRate);
                    r.setRw(rightEnd, mutationRate);
                } else {
                    r.setFw(rightEnd, mutationRate);
                    r.setRw(leftEnd, mutationRate);
                }

                // plot 2:
                numberOfMutations.add(r.getMutationsfw().size() + r.getMutationsrw().size());

                // Plot 3:
                numAllReads += 2;

                // fw
                if(fwGenomicRegions.size() == 1){ // read is NOT split
                    numNONSplit++;
                    if(r.getMutationsfw().isEmpty()) numNONSplitNOMm++;
                } else {
                    numSplit++;
                    if(r.getMutationsfw().isEmpty()) {
                        numSplitNOMm++;
                        // check fw regions; rw regions
                        for(Interval interval: fwGenomicRegions){
                            if(interval.getLength() < 5) break;
                        }
                        numSplitNOMsRegions++;
                    }

                }

                // rw
                if(rwGenomicRegions.size() == 1){
                    numNONSplit++;
                    if(r.getMutationsrw().isEmpty()) numNONSplitNOMm++;
                } else {
                    numSplit++;
                    if(r.getMutationsrw().isEmpty()) {
                        numSplitNOMm++;
                        for(Interval interval: rwGenomicRegions){
                            if(interval.getLength() < 5) break;
                        }
                        numSplitNOMsRegions++;
                    }
                }

                // New mutations distribution
                this.mutationPosition.addAll(r.getMutationsfw());
                this.mutationPosition.addAll(r.getMutationsrw());
            }
        }
    }

    /** get genomic regions for reads*/
    private ArrayList<Interval> regionVectors(int readStart, int readEnd, PreTranscript preTr){
        ArrayList<Interval> fwGenomicRegions = new ArrayList<>();

        // for(Exon ad: preTr.getExons()) System.out.println(Utils.GREEN + ad + Utils.RESET);

        for(Exon e: preTr.getExons()){
            // to locate the start (in which exon)
            if(readStart < e.getEndRelative()){ // special cases, 1. read == exon; 2. whole read is in one exon
                if(readEnd <= e.getEndRelative()){  // read completely in Exon
                    fwGenomicRegions.add(new Interval(e.getStartGenomic() + (readStart - e.getStartRelative()), e.getStartGenomic() + (readStart - e.getStartRelative()) + readLength));
                    break;
                }

                // otherwise read continues in the next exon
                fwGenomicRegions.add(new Interval(e.getStartGenomic() + (readStart-e.getStartRelative()), e.getEndGenomic()));

                e = preTr.getExons().higher(e);
                if(e != null){
                    while(readEnd > e.getEndRelative()){
                        fwGenomicRegions.add(new Interval(e.getStartGenomic(), e.getEndGenomic()));
                        e = preTr.getExons().higher(e);
                        if (e == null) System.err.println("End of read not found!");
                    }
                    fwGenomicRegions.add(new Interval(e.getStartGenomic(), e.getStartGenomic() + (readEnd-e.getStartRelative())));
                    break;
                } else System.err.println("End of read not found!");
            }
        }
        return fwGenomicRegions;
    }


    protected void createPlotsFrLength(){
        // Convert ArrayList to a comma-separated string
        StringBuilder data = new StringBuilder();
        for (long value : this.fragmentLengths) {
            data.append(value).append(",");
        }
        // Remove the last comma
        if (!data.isEmpty()) {
            data.setLength(data.length() - 1);
        }

        // Save data in a file
        try (FileWriter writer = new FileWriter("src/ReadSimulator/output/frLengths.csv")) {
            writer.write(data + "");
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            // Call the Python script and pass the data
            ProcessBuilder pb = new ProcessBuilder("python", "src/ReadSimulator/simulation_frlength.py");
            pb.inheritIO();  // Allows to see any output or errors from the Python script
            Process process = pb.start();
            int exitCode = process.waitFor();  // Wait for the script to finish
            if (exitCode == 0) {
                System.out.println("Plot generated successfully.");
            } else {
                System.out.println("Error: Python script did not run successfully.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    protected void createPlotsMutations(){
        // Convert ArrayList to a comma-separated string
        StringBuilder data = new StringBuilder();
        for (int value : this.numberOfMutations) {
            data.append(value).append(",");
        }
        // Remove the last comma
        if (!data.isEmpty()) {
            data.setLength(data.length() - 1);
        }

        // Save data in a file
        try (FileWriter writer = new FileWriter("src/ReadSimulator/output/mutations.csv")) {
            writer.write(data + "");
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            // Call the Python script and pass the data
            ProcessBuilder pb = new ProcessBuilder("python", "src/ReadSimulator/simulation_mutations.py");
            pb.inheritIO();  // This will allow you to see any output or errors from the Python script
            Process process = pb.start();
            int exitCode = process.waitFor();  // Wait for the script to finish
            if (exitCode == 0) {
                System.out.println("Plot generated successfully.");
            } else {
                System.out.println("Error: Python script did not run successfully.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    protected void createPlotsReadStatistics(){
        String[] values = new String[6];
        values[0] = numAllReads + "";
        values[1] = numNONSplit + "";
        values[2] = numNONSplitNOMm + "";
        values[3] = numSplit + "";
        values[4] = numSplitNOMm + "";
        values[5] = numSplitNOMsRegions + "";

        // Build the command to execute the Python script with arguments
        String[] command = new String[values.length + 2];
        command[0] = "python"; // Or "python3" depending on your environment
        command[1] = "src/ReadSimulator/simulation_read_statistics.py";
        System.arraycopy(values, 0, command, 2, values.length);

        System.out.println("Check python command");
        for(String s: command) System.out.println(s);

        try {
            // Execute the command
            ProcessBuilder pb = new ProcessBuilder(command);
            Process process = pb.start();

            // Capture the output from the Python script
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String line;
            while ((line = reader.readLine()) != null) {
                System.out.println(line);
            }

            // Wait for the process to complete and get the exit value
            int exitCode = process.waitFor();
            if (exitCode == 0) {
                System.out.println("Python script executed successfully.");
            } else {
                System.out.println("Python script encountered an error.");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    protected void createPlotsMutationPosition(){
        // Convert ArrayList to a comma-separated string
        StringBuilder data = new StringBuilder();
        for (long value : this.mutationPosition) {
            data.append(value).append(",");
        }
        // Remove the last comma
        if (!data.isEmpty()) {
            data.setLength(data.length() - 1);
        }

        // Save data in a file
        try (FileWriter writer = new FileWriter("src/ReadSimulator/output/mutationPosition.csv")) {
            writer.write(data + "");
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            // Call the Python script and pass the data
            ProcessBuilder pb = new ProcessBuilder("python", "src/ReadSimulator/simulation_mutation_position.py");
            pb.inheritIO();  // This will allow you to see any output or errors from the Python script
            Process process = pb.start();
            int exitCode = process.waitFor();  // Wait for the script to finish
            if (exitCode == 0) {
                System.out.println("Plot generated successfully.");
            } else {
                System.out.println("Error: Python script did not run successfully.");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}


