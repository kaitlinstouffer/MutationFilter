package com.example.mutationfilter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Dictionary;
import java.util.TreeMap;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/4/13
 * Time: 7:59 PM
 * To change this template use File | Settings | File Templates.
 *
 * Data Structure to Organize Files into Familial Relations
 */
public class FamilyDataGroup {

    //Instance Variables
    public String[] individuals;
    boolean[] chrReferences; // each reference is true if file lists locations with "chr"
    public ArrayList<ArrayList<Mutation>> possibleLocs;
    public ArrayList<Dictionary<String,ArrayList<Mutation>>> geneBanks;
    public ArrayList<Mutation> commonLocs;
    public TreeMap<String, ArrayList<Mutation>> commonGenes;
    // TO DO: implement flexible versions of these data structures

    // Constructors
    public FamilyDataGroup(String[] inputFileNames) {
        individuals = inputFileNames;
    }

    // Methods involving Family Data Groups and Instance Variables
    public static ArrayList<FamilyDataGroup> makeDataGroups(String inputFile) throws IOException {
        ArrayList<FamilyDataGroup> families = new ArrayList<FamilyDataGroup>();

        // Read in each line as family of tab separated individuals
        // Add set of Filenames to bank of families
        TabFileReader reader = new TabFileReader(inputFile);
        String[] files = reader.readColumnsArray();
        while (files != null) {
            FamilyDataGroup fdg = new FamilyDataGroup(files);
            families.add(fdg);
            files = reader.readColumnsArray();
        }

        return families;
    }

    /* Returns TreeMap of individuals, each with list of locations (chr:loc string)
     * that should be read in exon Read program.
     */
    public TreeMap<Integer, ArrayList<String>> getLocsToRead(int numOfIndividuals) {
        if (numOfIndividuals >= individuals.length) {
            return null;
        }

        // Examine possible Locs
        TreeMap<String, ArrayList<Integer>> preliminaryParse = new TreeMap<String,ArrayList<Integer>>();
        int numInd = 0;
        for (ArrayList<Mutation> aM : possibleLocs) {
            for (Mutation m : aM) {
                String s = m.chromLoc;
                if (chrReferences[numInd]) {
                    s = "chr" + m.chromLoc;
                }
                if (preliminaryParse.containsKey(s)) {
                    preliminaryParse.get(s).add(numInd);
                }
                else {
                    preliminaryParse.put(s,new ArrayList<Integer>());
                    preliminaryParse.get(s).add(numInd);
                }
            }
            numInd++;
        }

        // Make new TreeMap with mutations to read per individual
        TreeMap<Integer, ArrayList<String>> mutsToRead = new TreeMap<Integer,ArrayList<String>>();
        for (String s : preliminaryParse.keySet()) {
            int size = preliminaryParse.get(s).size();
            // Only consider the mutations that have enough, but not all individuals having them
            if ((size >= numOfIndividuals) && (size < individuals.length)) {
                ArrayList<Integer> withMut = preliminaryParse.get(s);
                for (int i = 0; i < individuals.length; i++) {
                    if (!withMut.contains(i)) {
                        if (mutsToRead.containsKey(i)) {
                            mutsToRead.get(i).add(s);
                        }
                        else {
                            mutsToRead.put(i, new ArrayList<String>());
                            mutsToRead.get(i).add(s);
                        }
                    }
                }
            }
        }
        return mutsToRead;
    }

    /* Sets commonLocs: an ArrayList of Mutations that are contained in all
     * individuals' possible locs lists as equal or unsequenced locations.
     * Mutations are equal if have same transcript, location, and base change.
     */
    public void setCommonLocs() {
        // keep TreeMap of all of locs and how many individuals possess them
        // note: only mutation objects that are sequenced are stored in commonLocs
        TreeMap<Mutation, Integer> commonLocTree = new TreeMap<Mutation,Integer>();
        ArrayList<Mutation> unsequencedMuts = new ArrayList<Mutation>();
        for (ArrayList<Mutation> am : possibleLocs) {
            for (Mutation m : am) {
                if (!m.isSequenced) {
                    unsequencedMuts.add(m);
                }
                else if (commonLocTree.containsKey(m)) {
                    commonLocTree.put(m,commonLocTree.get(m) + 1);
                }
                else {
                    commonLocTree.put(m,1);
                }
            }
        }
        // update counts in commonLocTree with unsequenced mutations
        // If more than one mutation exists at unsequenced location, all counts are updated
        for (Mutation m : unsequencedMuts) {
            for (Mutation mTree : commonLocTree.keySet()) {
                if (mTree.equalsOrUnsequenced(m)) {
                    commonLocTree.put(mTree,commonLocTree.get(mTree) + 1);
                }
            }
        }

        // Set Common Locs with all mutation that have count = # individuals
        // Note: Common locs may contain same mutations with different transcripts (8/1/2014)
        commonLocs = new ArrayList<Mutation> ();
        for (Mutation m : commonLocTree.keySet()) {
            if (commonLocTree.get(m).equals(individuals.length)) {
                commonLocs.add(m);
            }
        }

    }

    /* Sets commonGenes: if disease filtering is for recessive, non-consanguineous, then only
     * genes with at least two distinct mutations are retained.
     */
    public void setCommonGenes(Flags flags) {
        commonGenes = new TreeMap<String, ArrayList<Mutation>>();
        for (Mutation m : commonLocs) {
            if (!commonGenes.containsKey(m.gene)) {
                commonGenes.put(m.gene,new ArrayList<Mutation>());
            }
            commonGenes.get(m.gene).add(m);
        }

        // if recessive, non-consanguineous, then cut out all non-X genes that have < 2 mutations
        if (flags.RECESSIVE && (flags.CONSANGUINEOUS < 0)) {

            ArrayList<String> keysToKeep = new ArrayList<String>();
            for (String gene : commonGenes.keySet()) {
                if (commonGenes.get(gene).get(0).chromosome.equalsIgnoreCase("x")) {
                    keysToKeep.add(gene);
                }
                else if (commonGenes.get(gene).size() > 1) {
                    keysToKeep.add(gene);
                }
            }
            // Make new Dictionary
            TreeMap<String, ArrayList<Mutation>> newCommonGenes = new TreeMap<String, ArrayList<Mutation>>();
            for (String gene: keysToKeep) {
                newCommonGenes.put(gene, commonGenes.get(gene));
            }
            commonGenes = newCommonGenes;
        }
    }
}
