package org.molgenis;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class Chaperones {

    public static void main(String args[]) throws Exception {
        System.out.println("Starting...");

        // input files
        File chaperonesUniprotGenename = new File("/Users/joeri/git/dave/data/RAW/chaperones-uniprot-genename.txt");
        File bioGRID = new File("/Users/joeri/git/dave/data/RAW/BIOGRID-ORGANISM-Homo_sapiens-4.4.234.tab3.txt");

        // output files
        File chapIntBGLines = new File("/Users/joeri/git/dave/data/chap-biogrid.txt");
        File uniqGeneSymbols = new File("/Users/joeri/git/dave/data/interacting-with-chaperones-genenames.txt");
        File uniqUniprots = new File("/Users/joeri/git/dave/data/interacting-with-chaperones-uniprot.txt");

        // data structures
        Set<String> chapGeneSymbols = new TreeSet<>();
        Set<String> chapGeneUniprotIDs = new TreeSet<>();
        Set<String> chapGeneEnsGeneIDs = new TreeSet<>();

        Set<String> chapInteractorGeneSymbols = new TreeSet<>();
        Set<String> chapInteractorUniprotIDs = new TreeSet<>();
        Set<String> anyMatchingLines = new TreeSet<>();


        Scanner sc = new Scanner(chaperonesUniprotGenename);
        sc.nextLine(); //skip header
        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            String[] split = line.split("\t", -1);
            if (split.length != 3) {
                throw new Exception("Expected split 3 in line " + line);
            }
            String geneSymbol = split[0];
            String uniprot = split[1];
            String ensGeneID = split[2];

            if (chapGeneSymbols.contains(geneSymbol)) {
                System.out.println("dup: " + geneSymbol);
            }
            if (chapGeneUniprotIDs.contains(uniprot)) {
                System.out.println("dup: " + uniprot);
            }
            if (chapGeneEnsGeneIDs.contains(ensGeneID)) {
                System.out.println("dup: " + ensGeneID);
            }

            chapGeneSymbols.add(geneSymbol);
            chapGeneUniprotIDs.add(uniprot);
            chapGeneEnsGeneIDs.add(ensGeneID);
        }

        System.out.println("found " + chapGeneSymbols.size() + " unique gene symbols");
        System.out.println("found " + chapGeneUniprotIDs.size() + " unique uniprotIDs symbols");
        System.out.println("found " + chapGeneEnsGeneIDs.size() + " unique ensembl gene ids");


        /*
        indices:
        0:  #BioGRID Interaction ID
        1:  Entrez Gene Interactor A
        2:  Entrez Gene Interactor B
        3:  BioGRID ID Interactor A
        4:  BioGRID ID Interactor B
        5:  Systematic Name Interactor A
        6:  Systematic Name Interactor B
        7:  Official Symbol Interactor A
        8:  Official Symbol Interactor B
        9:  Synonyms Interactor A
        10: Synonyms Interactor B
        11: Experimental System
        12: Experimental System Type
        13: Author
        14: Publication Source
        15: Organism ID Interactor A
        16: Organism ID Interactor B
        17: Throughput
        18: Score
        19: Modification
        20: Qualifications
        21: Tags
        22: Source Database
        23: SWISS-PROT Accessions Interactor A
        24: TREMBL Accessions Interactor A
        25: REFSEQ Accessions Interactor A
        26: SWISS-PROT Accessions Interactor B
        27: TREMBL Accessions Interactor B
        28: REFSEQ Accessions Interactor B
        29: Ontology Term IDs
        30: Ontology Term Names
        31: Ontology Term Categories
        32: Ontology Term Qualifier IDs
        33: Ontology Term Qualifier Names
        34: Ontology Term Types
        35: Organism Name Interactor A
        36: Organism Name Interactor B
        */
        Scanner bgSc = new Scanner(bioGRID);
        bgSc.nextLine(); //skip header
        while (bgSc.hasNextLine()) {
            String line = bgSc.nextLine();
            String[] split = line.split("\t", -1);
            if (split.length != 37) {
                throw new Exception("Expected split 38 but was " + split.length + " in line:\n" + line);
            }


            String Official_Symbol_Interactor_A = split[7];
            String Official_Symbol_Interactor_B = split[8];
            String SWISS_PROT_Accessions_Interactor_A =  split[23];
            String SWISS_PROT_Accessions_Interactor_B =  split[26];


            if(chapGeneUniprotIDs.contains(Official_Symbol_Interactor_A) || chapGeneUniprotIDs.contains(SWISS_PROT_Accessions_Interactor_A))
            {
                chapInteractorGeneSymbols.add(Official_Symbol_Interactor_B);
                // in a few occasions, contain multiple IDs, e.g. B3EWG5|B3EWG6|B3EWG4|B3EWG3
                chapInteractorUniprotIDs.addAll(Arrays.asList(SWISS_PROT_Accessions_Interactor_B.split("\\|")));
                anyMatchingLines.add(line);
            }
            else if(chapGeneUniprotIDs.contains(Official_Symbol_Interactor_B) || chapGeneUniprotIDs.contains(SWISS_PROT_Accessions_Interactor_B))
            {
                chapInteractorGeneSymbols.add(Official_Symbol_Interactor_A);
                // in a few occasions, contain multiple IDs, e.g. B3EWG5|B3EWG6|B3EWG4|B3EWG3
                chapInteractorUniprotIDs.addAll(Arrays.asList(SWISS_PROT_Accessions_Interactor_A.split("\\|")));
                anyMatchingLines.add(line);
            }
        }

        System.out.println("unique: chaperone interactors by gene symbol: " + chapInteractorGeneSymbols.size());
        System.out.println("unique: chaperone interactors by uniprot ID: " + chapInteractorUniprotIDs.size());

        // write raw data for reference
        PrintWriter writer = new PrintWriter(chapIntBGLines, "UTF-8");
        for (String line : anyMatchingLines) {
            writer.write(line + System.lineSeparator());
        }
        writer.flush();
        writer.close();

        // write gene symbols
        PrintWriter writer2 = new PrintWriter(uniqGeneSymbols, "UTF-8");
        for (String gene : chapInteractorGeneSymbols) {
            writer2.write(gene + System.lineSeparator());
        }
        writer2.flush();
        writer2.close();

        // write uniprot IDs
        PrintWriter writer3 = new PrintWriter(uniqUniprots, "UTF-8");
        for (String gene : chapInteractorUniprotIDs) {
            writer3.write(gene + System.lineSeparator());
        }
        writer3.flush();
        writer3.close();

        System.out.println("...done!");

    }
}
