package ReadSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class FASTAIndex {
    private final HashMap<String, Long> chr2chrLength;
    private final HashMap<String, Long> chr2startIndexFASTA;
    private final ArrayList<String> chromosomes;
    private int lineLength;
    private int lineLengthInclNewline;

    public FASTAIndex(String findxFilePath) {
        this.chromosomes = new ArrayList<>();
        this.chr2startIndexFASTA = new HashMap<>();
        this.chr2chrLength = new HashMap<>();

        int count = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(findxFilePath))){
            String line;
            while ((line = br.readLine()) != null){
                String[] values = line.split("\t");
                if(count == 0) {
                    this.lineLength = Integer.parseInt(values[3]);
                    this.lineLengthInclNewline = Integer.parseInt(values[4]);
                }
                this.chromosomes.add(values[0]);
                this.chr2chrLength.put(values[0], Long.parseLong(values[1]));
                this.chr2startIndexFASTA.put(values[0], Long.parseLong(values[2]));

                count++;
            }
        }catch (IOException e) {
            e.printStackTrace();
        }
    }

    // Getter and Setter

    public ArrayList<String> getChromosomes() {
        return chromosomes;
    }

    public HashMap<String, Long> getChr2chrLength() {
        return chr2chrLength;
    }

    public HashMap<String, Long> getChr2startIndexFASTA() {
        return chr2startIndexFASTA;
    }

    public int getLineLength() {
        return lineLength;
    }

    public int getLineLengthInclNewline() {
        return lineLengthInclNewline;
    }
}
