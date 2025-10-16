package org.molgenis;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class VKGLProcessIntoProteinChanges {

    public static void main(String args[]) throws Exception {
        System.out.println("Starting...");

        // Input files
        // ** FOR ALL GENES, NEED TO RUN THREE TIMES, FOR EACH OF THESE FILES: **
        // --> protein-atlas-secreted-geneIDs-mane-uniprot (+.txt / +-withvariants.txt)
        // --> protein-atlas-intracellular-geneIDs-mane-uniprot-random2000 (+.txt / +-withvariants.txt)
        // --> protein-atlas-membrane-geneIDs-mane-uniprot-random2000 (+.txt / +-withvariants.txt)
        String inputGenes = "protein-atlas-membrane-geneIDs-mane-uniprot-random2000";
        File proteinSetWithMane = new File("/Users/joeri/git/dave1/data/" + inputGenes + ".txt");
        String vkglMissenseFileName = "VKGL_apr2024_annot_missense_nodup_b38.vcf";
        File vkglMissenseLocation = new File("/Users/joeri/git/dave1/data/" + vkglMissenseFileName + ".zip");

        // Output files
        // Directory to write all separate variants files to, creating 1 additional directory per gene:
        File outputBaseDir = new File("/Users/joeri/git/dave1/data/genes/");
        // File with genes having 1+ variants
        File genesWithVariants = new File("/Users/joeri/git/dave1/data/" + inputGenes + "-withvariants.txt");

        System.out.println("Loading gene/transcript file...");
        HashMap<String, String> geneNameToTranscript = new HashMap<>();
        HashMap<String, String> geneNameToUniprot = new HashMap<>();
        Set<String> geneNamesWithVariants = new TreeSet<>();
        Set<String> geneStableIdUniqueCheck = new TreeSet<>();
        Set<String> transcriptUniqueCheck = new TreeSet<>();
        Scanner sc = new Scanner(proteinSetWithMane);
        System.out.println("Skipping header: " + sc.nextLine());
        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            String[] split = line.split("\t", -1);
            String geneStableID = split[0];
            String transcript = split[1];
            String uniprot = split[2];
            String geneName = split[3];
            if (split.length != 4) {
                throw new Exception("Expected split 3 in line " + line);
            }
            if (geneNameToTranscript.containsKey(geneName)) {
                throw new Exception("Gene name duplicate " + geneName);
            }
            geneNameToTranscript.put(geneName, transcript);
            geneNameToUniprot.put(geneName, uniprot);
            if(geneStableIdUniqueCheck.contains(geneStableID))
            {
                throw new Exception("Duplicate gene stable ID: " + geneStableID);
            }
            geneStableIdUniqueCheck.add(geneStableID);
            if(transcriptUniqueCheck.contains(transcript))
            {
                throw new Exception("Duplicate transcript stable ID: " + transcript);
            }
            transcriptUniqueCheck.add(transcript);
        }
        sc.close();
        System.out.println("...done");

        System.out.println("Pre-loading all lines from ZIP file...");
        Map<String, List<String>> trancriptToLines = new HashMap<>();
        ZipInputStream zipStream = new ZipInputStream(new FileInputStream(vkglMissenseLocation));
        ZipEntry entry;
        while ((entry = zipStream.getNextEntry()) != null) {
            System.out.println("ZIP entry name: " + entry.getName());
            if (!entry.getName().equals(vkglMissenseFileName)) {
                continue;
            } else {
                System.out.println("ZIP entry for " + vkglMissenseFileName + " found, processing...");
            }
            Scanner s = new Scanner(zipStream);
            String line;
            while (s.hasNextLine()) {
                line = s.nextLine();

                String[] lineSplit = line.split("\t", -1);
                if (lineSplit.length != 8) {
                    throw new Exception("line split expecting 8 elements for " + line);
                }

                String[] infoSplit = lineSplit[7].split(";", -1);
                if (infoSplit.length != 2) {
                    throw new Exception("info split expecting 2 elements for " + line);
                }

                String[] csqSplit = infoSplit[1].split(",", -1);
                if (csqSplit.length == 0) {
                    throw new Exception("csq split expecting 1 or more elements for " + line);
                }

                for (String csq : csqSplit) {
                    String[] csqParts = csq.split("\\|", -1);
                    // CSQ=A|missense_variant|MODERATE|AGRN|ENSG00000188157|Transcript|ENST00000379370|protein_coding|22/36||||3783|3733|1245|V/M|Gtg/Atg|||1||HGNC|329
                    if (csqParts.length != 23) {
                        throw new Exception("csq parts split expecting 23 elements but was " + csqParts.length + " for " + line);
                    }

                    // MANE transcript
                    String transcript = csqParts[6];

                    if (trancriptToLines.containsKey(transcript)) {
                        trancriptToLines.get(transcript).add(line);
                    } else {
                        List<String> lines = new ArrayList<>();
                        lines.add(line);
                        trancriptToLines.put(transcript, lines);
                    }
                }
            }
        }
        zipStream.close();
        System.out.println("...done");

        System.out.println("Iterating over all genes in ZIP contents...");
        for (String geneName : geneNameToTranscript.keySet()) {
            System.out.println(" gene: " + geneName);

            Map<String, String> protChangesToClsf = new TreeMap<>();
            Map<String, String> protChangesToChromPosRefAlt = new TreeMap<>();
            Map<String, String> protChangesToAALoc = new TreeMap<>();

            List<String> lines = trancriptToLines.get(geneNameToTranscript.get(geneName));
            if (lines == null) {
                System.out.println("No lines for " + geneName + "!");
                continue;
            }
            for (String line : lines) {

                String[] lineSplit = line.split("\t", -1);
                String[] infoSplit = lineSplit[7].split(";", -1);
                String[] csqSplit = infoSplit[1].split(",", -1);

                for (String csq : csqSplit) {
                    String[] csqParts = csq.split("\\|", -1);

                    // sanity-check: exact match to MANE transcript
                    if (!csqParts[6].equals(geneNameToTranscript.get(geneName))) {
                        continue;
                    }

                    String proteinPos = csqParts[14];

                    // some variants are interpreted in a non-coding context in which case proteinPos is empty
                    if (proteinPos.isBlank()) {
                        continue;
                    }
                    String[] proteinChangeSplit = csqParts[15].split("/");
                    // synonymous variants have only 1 letter ('A') instead of 'A/B'
                    if (proteinChangeSplit.length != 2) {
                        continue;
                    }

                    String geneSymbol = csqParts[3];
                    // as per FoldX notation: WT residue, chain, residue number, mutant residue (e.g. "CA1490Y;")
                    String proteinChange = proteinChangeSplit[0] + "A" + proteinPos + proteinChangeSplit[1];
                    if (proteinChange.contains("=")) {
                        System.out.println("Synonymous variant, skipping " + proteinChange);
                        continue;
                    }
                    if (proteinChange.contains("*")) {
                        System.out.println("Stop gain, skipping " + proteinChange);
                        continue;
                    }
                    if(proteinChange.contains("-"))
                    {
                        System.out.println("Double notation (e.g. NVA287-288NI), skipping " + proteinChange);
                        continue;
                    }
                    String geneProt = geneSymbol + "\t" + proteinChange;
                    String clf = infoSplit[0].replace("VKGL=", "");

                    // unless there is a different protein notation for a transcript, the notation is duplicate
                    // filter here by checking if gene-prot combination is already present
                    // at the same time check for any potential conflicting interpretations for this notation
                    if (protChangesToClsf.containsKey(geneProt)) {
                        String prevClf = protChangesToClsf.get(geneProt);
                        if (!prevClf.equals(clf)) {
                            System.out.println("Conflicting classification for key '" + geneProt + "': " + prevClf + " vs " + clf + " at " + line);
                            protChangesToClsf.put(geneProt, "CF");
                        }
                    } else {
                        // first time add, also store genome coordinates
                        protChangesToClsf.put(geneProt, clf);
                        protChangesToChromPosRefAlt.put(geneProt, (lineSplit[0].replace("chr", "") + "\t" + lineSplit[1] + "\t" + lineSplit[3] + "\t" + lineSplit[4]));
                        protChangesToAALoc.put(geneProt, proteinPos);
                    }
                }
            }


        /*
        Write results
         */
            if (protChangesToClsf.isEmpty()) {
                continue;
            }
            File geneDir = new File(outputBaseDir + "/" + geneName);
            if (!geneDir.exists()) {
                geneDir.mkdirs();
            }
            PrintWriter writer = new PrintWriter(outputBaseDir + "/" + geneName + "/VKGL_apr2024_protForFolding.tsv", "UTF-8");
            writer.write("Assembly" + "\t" + "Chrom" + "\t" + "Pos" + "\t" + "Ref" + "\t" + "Alt" + "\t" + "Gene" + "\t" + "ProtChange" + "\t" + "AALoc" + "\t" + "Classification" + System.lineSeparator());
            for (String key : protChangesToClsf.keySet()) {
                writer.write("GRCh38" + "\t" + protChangesToChromPosRefAlt.get(key) + "\t" + key + "\t" + protChangesToAALoc.get(key) + "\t" + protChangesToClsf.get(key) + System.lineSeparator());
            }
            writer.flush();
            writer.close();
            geneNamesWithVariants.add(geneName);
        }
        System.out.println("...done");


        System.out.println("Writing file which genes had >1 variants...");
        PrintWriter writer2 = new PrintWriter(genesWithVariants, "UTF-8");
        writer2.write("Transcript stable ID\tGene name\tUniProtKB/Swiss-Prot ID" + System.lineSeparator());
        for (String geneName : geneNamesWithVariants) {
            writer2.write(geneNameToTranscript.get(geneName) + "\t" + geneName + "\t" + geneNameToUniprot.get(geneName) + System.lineSeparator());
        }
        writer2.flush();
        writer2.close();
        System.out.println("...done");
    }
}
