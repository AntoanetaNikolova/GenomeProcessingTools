package ExonSkipping;

import BaseComponents.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Store all relevant info from the GTF file
 * @Input: a GTF filename
 * @Output: Table with all ES-SE events - in the rows; wirth 14 columns with info
 */

public class ExonSkipping {
    private Genome genome;

    /**
     * @param inputFilePath reads file, iterates it and saves the described gene structures
     * 1. chromosome
     * 2. source (gene biotype)
     * 3. type of feature (gene, transcript, CDS...)
     * 4. start
     * 5. end
     * 6. score *irrelevant*
     * 7. strand +/-
     * 8. frame *irrelevant*
     * 9. attribute - their order and number of elements *can* change
     */
    public ExonSkipping(String inputFilePath) {
        genome = new Genome();

        try (BufferedReader br = new BufferedReader(new FileReader(inputFilePath))) {
            String line;
            while ((line = br.readLine()) != null) {

                int firstTab = line.indexOf("\t");
                if (firstTab != -1) {
                    int secondTab = line.indexOf("\t", firstTab + 1);
                    if (secondTab != -1) {
                        String afterSecondTab = line.substring(secondTab + 1);
                        if (afterSecondTab.startsWith("CDS")){
                            String[] values = line.split("\t");
                            String[] attributes = values[8].split("; ");
                            CDS newCDS = new CDS(Integer.parseInt(values[3]), Integer.parseInt(values[4]));

                            Map<String, String> mappedAttributes = Arrays.stream(attributes)
                                    .filter(e -> e.startsWith("gene") || e.startsWith("transcript") || e.startsWith("protein") || e.startsWith("ccds"))  // filter out only attributes we are interested in
                                    .map(attr -> attr.split(" ", 2))
                                    .filter(pair -> pair.length == 2)
                                    .collect(Collectors.toMap(pair -> pair[0].trim(), pair -> pair[1].replace("\"", "").replace(";", ""), (existing, replacement) -> replacement));


                            String currentGeneID = mappedAttributes.getOrDefault("gene_id", "");
                            if(currentGeneID.isEmpty()){
                                currentGeneID = attributes[0].split(" ")[2].replace("\"", "").replace(";", "");
                            }
                            String currentGeneName = mappedAttributes.getOrDefault("gene_name", "");
                            String currentTranscriptName = mappedAttributes.getOrDefault("transcript_name", "");
                            String currentTranscriptID = mappedAttributes.getOrDefault("transcript_id", "");
                            String currentCDSID = mappedAttributes.getOrDefault("protein_id", mappedAttributes.getOrDefault("ccds_id", ""));
                            newCDS.setID(currentCDSID);

                            if (!genome.getId2Gene().containsKey(currentGeneID)) {  // If gene origin of CDS not in genome
                                // Define new Gene and add it; If Gene new, then also transcript new
                                Gene geneToAdd = new Gene(values[0], values[6]);
                                geneToAdd.setID(currentGeneID);
                                geneToAdd.setName(currentGeneName);

                                Transcript newTranscript = new Transcript(currentTranscriptID, currentTranscriptName);
                                newTranscript.addCDS(newCDS);

                                // Transcript to gene
                                geneToAdd.addTranscript(newTranscript);

                                // Add the new gene to genome.
                                genome.getId2Gene().put(currentGeneID, geneToAdd);
                                genome.getId2transcript().put(currentTranscriptID, newTranscript);
                                genome.addGene(geneToAdd);
                            } else {  // Gene is already in genome
                                Gene geneOrigin = genome.getId2Gene().get(currentGeneID);
                                Utils.checkIfNull(geneOrigin, new Gene());

                                // Check existence of transcript:
                                Transcript transcriptToUpdate;
                                if (!genome.getId2transcript().containsKey(currentTranscriptID)) {
                                    transcriptToUpdate = new Transcript(currentTranscriptID, currentTranscriptName);
                                    transcriptToUpdate.addCDS(newCDS);
                                    genome.getId2transcript().put(currentTranscriptID, transcriptToUpdate);
                                    geneOrigin.addTranscript(transcriptToUpdate);
                                } else {
                                    // it exists => just add element
                                    transcriptToUpdate = genome.getId2transcript().get(currentTranscriptID);
                                    Utils.checkIfNull(transcriptToUpdate, new Transcript());

                                    transcriptToUpdate.addCDS(newCDS);
                                }
                            }
                        }
                    }
                }
            }


            for (Gene g : genome.getAllGenes()) {
                g.setN_ports(g.getAllTranscripts().size());
                g.setN_trans(g.getAllTranscripts().size());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * For each gene define all introns by iterating each transcript in the gene
     * @return the ES-SE string
     */
    public String defineAndOrganizeIntrons() {
        for (Gene gene : genome.getAllGenes()) {
            ArrayList<Intron> introns = new ArrayList<>();
            HashMap<Integer, ArrayList<Intron>> organizedIntrons = new HashMap<>();

            for (Transcript transcript: gene.getAllTranscripts()) {
                ArrayList<Integer> startingPositions = new ArrayList<>();

                CDS end5 = null;
                for (CDS end3 : transcript.getAllCDS()) {  // CDS entries are automatically sorted
                    if (end5 != null) {
                        Intron intron = new Intron(end5, end3);
                        intron.setOriginTranscript(transcript);
                        transcript.getAllIntrons().add(intron);
                        introns.add(intron);

                        int startIndex = intron.getStartGenomic();
                        startingPositions.add(startIndex);

                        if(!organizedIntrons.containsKey(startIndex)){
                            organizedIntrons.put(startIndex, new ArrayList<>());
                            organizedIntrons.get(startIndex).add(intron);
                        } else organizedIntrons.get(startIndex).add(intron);
                    }
                    end5 = end3;
                }

                Collections.sort(startingPositions);
                transcript.setIntronStartingPositions(startingPositions);
            }

            gene.setAllIntrons(introns);
            gene.setOrganizedIntrons(organizedIntrons);

        }

        return this.ES_SE();
    }


    /**
     * Uses the data structure with organized introns;
     * Has a list with all introns gro gene;
     * Iterate through all genes and build output
     */
    private String ES_SE() {
        genome.getAllGenes().forEach(gene -> {

            List<Intron> uncheckedIntrons = gene.getAllIntrons().stream()
                    .filter(intron -> !intron.getChecked())
                    .toList();

            uncheckedIntrons.forEach(intron -> processIntron(gene, intron));
        });

        String output = "id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n"
                + genome.getExonSpliceEvents().stream()
                .map(ES_SE::toString)
                .collect(Collectors.joining("\n"));

        return output;
    }

    private void processIntron(Gene gene, Intron intron) {
        if (intron.getChecked()) return;

        int SVIntervalStart = intron.getStartGenomic();
        int SVIntervalEnd = intron.getEndGenomic();

        ArrayList<Intron> SV = new ArrayList<>();
        ArrayList<Intron> WT = new ArrayList<>();

        ArrayList<Intron> potentialTranscripts = gene.getOrganizedIntrons().get(SVIntervalStart);

        if (potentialTranscripts.size() == 1) {
            potentialTranscripts.getFirst().setChecked();
            return;
        }

        for (Intron potentialEvent : potentialTranscripts) {
            if (potentialEvent.getEndGenomic() > SVIntervalEnd) continue;

            if (potentialEvent.getEndGenomic() == SVIntervalEnd) {
                potentialEvent.setChecked();
                SV.add(potentialEvent);
            } else {

                if(evaluateWT(gene, potentialEvent, SVIntervalEnd)){
                    WT.add(potentialEvent);
                }
            }
        }

        if (!WT.isEmpty()) {
            ArrayList<Intron> WT_introns = determineWTIntronsInSVInterval(WT, SVIntervalStart, SVIntervalEnd);
            recordExonSpliceEvent(gene, SV, WT, WT_introns, SVIntervalStart, SVIntervalEnd);
        }
    }

    private boolean evaluateWT(Gene gene, Intron potentialIntron, int SV_interval_end){
        Transcript transcript = potentialIntron.getOriginTranscript();
        String transcriptID = transcript.getID();

        int index = Collections.binarySearch(transcript.getIntronStartingPositions(), potentialIntron.getStartGenomic());
        if (index < 0) {
            System.err.println("STARTING INDEX NOT FOUND!");
            return false;
        }

        A: while (index < transcript.getIntronStartingPositions().size() && transcript.getIntronStartingPositions().get(index) < SV_interval_end) {
            int startPosition = transcript.getIntronStartingPositions().get(index);
            ArrayList<Intron> intronsAtPosition = gene.getOrganizedIntrons().get(startPosition);

            for (Intron candidate : intronsAtPosition) {
                if (!candidate.getOriginTranscript().getID().equals(transcriptID)) continue;

                // Candidate has the same transcript iD -> happens ones in a gene
                if (candidate.getEndGenomic() == SV_interval_end) {  // Defines event
                    return true;
                } else if (candidate.getEndGenomic() < SV_interval_end) {
                    index++;
                    continue A;
                } else
                    return false;
            }

            index++;
        }

        return false;
    }

    private void recordExonSpliceEvent(Gene gene, ArrayList<Intron> SV, ArrayList<Intron> WT, ArrayList<Intron> WT_introns, int SV_interval_start, int SV_interval_end) {
        ES_SE event = new ES_SE(gene.getID(), gene.getName(), gene.getChromosome(),
                gene.getStrandDirection().charAt(0), gene.getN_ports(),
                gene.getN_trans(), SV, WT, WT_introns);

        int minSkippedExons = calculateSkippedExons(WT, SV_interval_start, SV_interval_end, true);
        event.setMin_skipped_exon(minSkippedExons);
        int maxSkippedExons = calculateSkippedExons(WT, SV_interval_start, SV_interval_end, false);
        event.setMax_skipped_exon(maxSkippedExons);
        event.setMin_skipped_bases(calculateSkippedBases(WT, SV_interval_start, SV_interval_end, true));
        event.setMax_skipped_bases(calculateSkippedBases(WT, SV_interval_start, SV_interval_end, false));

        genome.getExonSpliceEvents().add(event);
    }

    private int calculateSkippedExons(ArrayList<Intron> WTs, int SV_startIndex, int SV_endIndex, boolean isMin) {
        return WTs.stream()
                .mapToInt(wt -> {
                    Transcript transcript = wt.getOriginTranscript();
                    return (int) transcript.getAllIntrons().stream()
                            .filter(intron -> intron.getStartGenomic() >= SV_startIndex && intron.getEndGenomic() <= SV_endIndex)
                            .count() - 1;
                })
                .reduce(isMin ? Integer.MAX_VALUE : Integer.MIN_VALUE, isMin ? Integer::min : Integer::max);
    }

    private int calculateSkippedBases(ArrayList<Intron> WTs, int SV_startIndex, int SV_endIndex, boolean isMin) {
        return WTs.stream()
                .mapToInt(wt -> {
                    Transcript transcript = wt.getOriginTranscript();
                    // Get filtered introns
                    List<Intron> filteredIntrons = transcript.getAllIntrons().stream()
                            .filter(intron -> intron.getStartGenomic() >= SV_startIndex && intron.getEndGenomic() < SV_endIndex)
                            .toList();

                    // Sum lengths of the CDS at the 3' end
                    int cdsLengthSum = filteredIntrons.stream()
                            .mapToInt(intron -> intron.getCds3End().getLength())
                            .sum();

                    // Count the number of filtered introns
                    int filteredIntronCount = filteredIntrons.size();

                    // Return the total sum of CDS length + filtered intron count
                    return cdsLengthSum + filteredIntronCount -1;
                })
                .reduce(isMin ? Integer.MAX_VALUE : Integer.MIN_VALUE, isMin ? Integer::min : Integer::max);
    }

    private ArrayList<Intron> determineWTIntronsInSVInterval(ArrayList<Intron> WTs, int SV_startIndex, int SV_endIndex) {
        Map<String, Intron> uniqueIntronsMap = new HashMap<>();

        WTs.forEach(wt -> {
            Transcript transcript = wt.getOriginTranscript();
            transcript.getAllIntrons().stream()
                    .filter(intron -> intron.getStartGenomic() >= SV_startIndex && intron.getEndGenomic() <= SV_endIndex)
                    .forEach(intron -> uniqueIntronsMap.put(
                            intron.getStartGenomic() + ":" + intron.getEndGenomic(),
                            intron
                    ));
        });

        // Return unique introns as a list
        return new ArrayList<>(uniqueIntronsMap.values());
    }


    // Getter and Setter
    public Genome getGenome() {
        return genome;
    }
}