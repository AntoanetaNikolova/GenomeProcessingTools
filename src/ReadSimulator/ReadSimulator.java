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
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.NormalDistribution;

public class ReadSimulator {
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

    private final String fwFastqPath;
    private final String rwFastqPath;
    private final String readMappingsInfoPath;


    public ReadSimulator(int length, double frlength, double SD, String readcounts , double mutationRate, String fasta, String fixd, String gtf, String outputPath) {
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

        this.fwFastqPath = outputPath + File.separator + "fw.fastq";
        this.rwFastqPath = outputPath  + File.separator  + "rw.fastq";
        this.readMappingsInfoPath = outputPath + File.separator + "read.mappinginfo";

        // add header
        try (BufferedWriter mappingsWriter = new BufferedWriter(new FileWriter(readMappingsInfoPath, true))){
            mappingsWriter.write("readid\tchr\tgene\ttranscript\tfw_regvec\trw_regvec\tt_fw_regvec\tt_rw_regvec\tfw_mut\trw_mut\n");
        }catch (IOException e) {
            System.err.println("Error writing in files: " + e.getMessage());
            e.printStackTrace();
        }
    }

    /** create a ReadCollection object for each line in readcounts-file */
    protected void defineReadCollections() {
        ReadCollection collection;
        try (BufferedReader br = new BufferedReader(new FileReader(readCountsFilePath))) {
            String line;
            while ((line = br.readLine()) != null){
                String[] elements = line.split("\t");
                if(elements[2].equals("count")) continue;
                if (elements.length != 3) System.err.println("Error! Line in readcounts-file with less/more than 3 elements");

                collection = new ReadCollection(elements[0], elements[1], Integer.parseInt(elements[2]));
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
        // declare variables to avoid excessive memory allocation
        int readCounter = 0;
        StringBuilder outputFwFASTQ;
        StringBuilder outputRwFASTQ;
        StringBuilder outputReadMapping;
        ArrayList<Interval> fwGenomicRegions;
        ArrayList<Interval> rwGenomicRegions;

        for(PreTranscript preTr: transcriptId2preTranscript.values()){
            // 1. define relative coordinates exons
            int transcriptLength = 0;
            int startRelativeIndex = 0;

            for(Exon sortedExon: preTr.getExons()){
                sortedExon.setStartRelative(startRelativeIndex);
                int endRelative = startRelativeIndex + sortedExon.getLength();
                sortedExon.setEndRelative(endRelative);
                startRelativeIndex = endRelative;
                transcriptLength += sortedExon.getLength();
            }

            // if the length of the transcript is smaller than the read length => not possible to extract reads
            preTr.setTranscriptLength(transcriptLength, readLength);
            if(preTr.getTooShort()) {
                transcriptID2ReadCollection.remove(preTr.getTranscriptID());  // there is no need to iterate it, no reads will be created from it
                continue;
            }

            // 2. extract pre-transcript sequence from FASTA using FASTA Index
            long preTranscriptStartGenomic = 0;
            long preTranscriptEndGenomic = 0;

            // 1-based and end-inclusive indexes
            preTranscriptStartGenomic = preTr.getExons().getFirst().getStartGenomic();
            preTranscriptEndGenomic = preTr.getExons().getLast().getEndGenomic() - 1;

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

            // 3. modify pre-transcript sequence - remove introns
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

            // 4. get the corresponding ReadCollection object; for each read-object define the fragment and read indexes
            ReadCollection readColCurTranscript = transcriptID2ReadCollection.get(preTr.getTranscriptID());
            Utils.checkIfNull(readColCurTranscript, new ReadCollection());
            if(readColCurTranscript == null) continue;

            readColCurTranscript.setChromosome(preTr.getChromosome());

            Read r;
            // Add
            outputFwFASTQ = new StringBuilder();
            outputRwFASTQ = new StringBuilder();
            outputReadMapping = new StringBuilder();

            for(int i = 0; i < readColCurTranscript.getNumReads(); i++){
                r = new Read(readLength);

                // determine fragment length
                long fragmentLength = Math.round(distribution.sample());
                while(fragmentLength < readLength || fragmentLength > preTr.getTranscriptLength()) {
                    fragmentLength = Math.round(distribution.sample());
                }

                // determine fragment start
                if(fragmentLength == preTr.getTranscriptLength()) {
                    r.setFragmentStartRelative(0, (int)fragmentLength, (int)preTr.getTranscriptLength(), preTr.getStrandDirection());
                } else {
                    r.setFragmentStartRelative(Math.max(0, new Random().nextInt((int)(preTr.getTranscriptLength() - fragmentLength))), (int)fragmentLength, (int)preTr.getTranscriptLength(), preTr.getStrandDirection());
                }

                // get genomic region vectors
                fwGenomicRegions = regionVectors(r.getFwStartRelative(), r.getFwEndRelative(), preTr);
                rwGenomicRegions = regionVectors(r.getRwStartRelative(), r.getRwEndRelative(), preTr);

                // set Fw and Rw
                String leftEnd = preTranscriptSequence.substring(r.getFwStartRelative(), r.getFwEndRelative());
                String rightEnd = Utils.getReverseComplement(new StringBuilder(preTranscriptSequence.substring(r.getRwStartRelative(), r.getRwEndRelative())));

                if(preTr.getStrandDirection().equals("+")){
                    r.setFw(leftEnd, mutationRate);
                    r.setRw(rightEnd, mutationRate);
                } else {
                    r.setFw(rightEnd, mutationRate);
                    r.setRw(leftEnd, mutationRate);
                }

                if(preTr.getStrandDirection().equals("-")){
                    outputFwFASTQ.append("@" + readCounter + "\n" + r.getFw() + "\n" + "+" + readCounter + "\n" + "I".repeat(r.getReadLength()) + "\n");
                    outputRwFASTQ.append("@" + readCounter + "\n" + r.getRw() + "\n" + "+" + readCounter + "\n" + "I".repeat(r.getReadLength()) + "\n");
                    outputReadMapping.append(readCounter + "\t" + readColCurTranscript.getChromosome() + "\t" + readColCurTranscript.getOriginGeneID() + "\t" + readColCurTranscript.getOriginTranscriptID() +
                            "\t" + rwGenomicRegions.stream().map(Interval::toString).collect(Collectors.joining("|")) + "\t" + fwGenomicRegions.stream().map(Interval::toString).collect(Collectors.joining("|")) + "\t" + r.getRwTrStartRelative() + "-" + r.getRwTrEndRelative() + "\t" + r.getFwTrStartRelative() + "-" + r.getFwTrEndRelative() + "\t" +
                            r.getMutationsfw().stream().map(String::valueOf).collect(Collectors.joining(", ")) + "\t" +
                            r.getMutationsrw().stream().map(String::valueOf).collect(Collectors.joining(", ")) + "\n");
                } else {
                    outputFwFASTQ.append("@" + readCounter + "\n" + r.getFw() + "\n" + "+" + readCounter + "\n" + "I".repeat(r.getReadLength()) + "\n");
                    outputRwFASTQ.append("@" + readCounter + "\n" + r.getRw() + "\n" + "+" + readCounter + "\n" + "I".repeat(r.getReadLength()) + "\n");
                    outputReadMapping.append(readCounter + "\t" + readColCurTranscript.getChromosome() + "\t" + readColCurTranscript.getOriginGeneID() + "\t" + readColCurTranscript.getOriginTranscriptID() +
                            "\t" + fwGenomicRegions.stream().map(Interval::toString).collect(Collectors.joining("|")) + "\t" + rwGenomicRegions.stream().map(Interval::toString).collect(Collectors.joining("|")) + "\t" + r.getFwTrStartRelative() + "-" + r.getFwTrEndRelative() + "\t" + r.getRwTrStartRelative() + "-" + r.getRwTrEndRelative() + "\t" +
                            r.getMutationsfw().stream().map(String::valueOf).collect(Collectors.joining(", ")) + "\t" +
                            r.getMutationsrw().stream().map(String::valueOf).collect(Collectors.joining(", ")) + "\n");
                }
                readCounter++;
            }

            // 5. output format
            addStringBuilders(outputFwFASTQ, outputRwFASTQ, outputReadMapping);
        }
    }

    /** append output string */
    private void addStringBuilders(StringBuilder fwQ, StringBuilder rwQ, StringBuilder readM){
        try (
                BufferedWriter fwWriter = new BufferedWriter(new FileWriter(fwFastqPath, true));
                BufferedWriter rwWriter = new BufferedWriter(new FileWriter(rwFastqPath, true));
                BufferedWriter mappingsWriter = new BufferedWriter(new FileWriter(readMappingsInfoPath, true))
        ){
            fwWriter.write(fwQ.toString());
            rwWriter.write(rwQ.toString());
            mappingsWriter.write(readM.toString());
        }catch (IOException e) {
            System.err.println("Error writing in files: " + e.getMessage());
            e.printStackTrace();
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
}

