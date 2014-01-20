package com.example.mutationfilter;

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
    public String[] bamIndividuals; // only used if bam files are specified
    boolean[] chrReferences; // each reference is true if file lists locations with "chr"
    public ArrayList<ArrayList<Mutation>> possibleLocs;
    public ArrayList<Mutation> commonLocs;
    public TreeMap<String, ArrayList<Mutation>> commonGenes;
    public boolean consanguineous;
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
        reader.close();

        return families;
    }

    // Sets bam file names
    public static void setBamFileNames(ArrayList<FamilyDataGroup> families, String inputFile) throws IOException {
        TabFileReader reader = new TabFileReader(inputFile);
        String[] files = reader.readColumnsArray();
        int count = 0;
        while (files != null) {
            families.get(count).bamIndividuals = files;
            files = reader.readColumnsArray();
            count++;
        }
        if (count != families.size()) {
            System.out.println("WARNING: not enough bam files specified in " + inputFile);
        }
        reader.close();
    }

    // Sets consanguinity flag
    public static void setConsanguinity(ArrayList<FamilyDataGroup> families, Flags flags) {
        for (FamilyDataGroup fdg : families) {
            fdg.consanguineous = false;
        }
        // if no consanguineous flag set, assume no families are consanguineous
        if (flags.CONSANGUINEOUS >= 0) {
            if (flags.CONSANGUINEOUS_FAMILIES != null) {
                for (Integer i : flags.CONSANGUINEOUS_FAMILIES) {
                    families.get(i).consanguineous = true;
                }
            }
            // apply consanguineous flag to all families if cFam was not specified
            else {
                for (FamilyDataGroup fdg : families) {
                    fdg.consanguineous = true;
                }
            }
        }
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
                            if (!mutsToRead.get(i).contains(s))
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
        TreeMap<Mutation, boolean[]> commonLocTree = new TreeMap<Mutation,boolean[]>();
        TreeMap<Mutation, ArrayList<Integer>> unsequencedMuts = new TreeMap<Mutation, ArrayList<Integer>>();
        int numIndividual = 0;
        for (ArrayList<Mutation> am : possibleLocs) {
            for (Mutation m : am) {
                if (!m.isSequenced) {
                    if (!unsequencedMuts.containsKey(m)) {
                        unsequencedMuts.put(m,new ArrayList<Integer>());
                    }
                    unsequencedMuts.get(m).add(numIndividual);
                }
                else if (commonLocTree.containsKey(m)) {
                    commonLocTree.get(m)[numIndividual] = true;
                }
                else {
                    commonLocTree.put(m,new boolean[individuals.length]);
                    commonLocTree.get(m)[numIndividual] = true;
                }
            }
            numIndividual++;
        }

        // update counts in commonLocTree with unsequenced mutations
        // If more than one mutation exists at unsequenced location, all counts are updated
        for (Mutation m : unsequencedMuts.keySet()) {
            for (Mutation mSeq : commonLocTree.keySet()) {
                if (m.equalsOrUnsequenced(mSeq)) {
                    for (Integer in : unsequencedMuts.get(m)) {
                        commonLocTree.get(mSeq)[in] = true;
                    }
                }
            }
        }

        // Set Common Locs with all mutation that have count = # individuals
        // Note: Common locs may contain same mutations with different transcripts (8/1/2014)
        commonLocs = new ArrayList<Mutation> ();
        for (Mutation m : commonLocTree.keySet()) {
            boolean[] temp = commonLocTree.get(m);
            boolean addToCommon = true;
            for (boolean b : temp) {
                if (!b) {
                    addToCommon = false;
                    break;
                }
            }
            if (addToCommon) {
                commonLocs.add(m);
            }
        }

        // TODO: REMOVE
        System.out.println("size of common locs is " + commonLocs.size());

    }

    /* Sets commonGenes: if disease filtering is for recessive, non-consanguineous, then only
     * genes with at least two distinct mutations are retained.
     */
    public void setCommonGenes(Flags flags) {
        commonGenes = new TreeMap<String, ArrayList<Mutation>>();
        TreeMap<String, ArrayList<Mutation>> mutCount = new TreeMap<String, ArrayList<Mutation>>();
        for (Mutation m : commonLocs) {
            if (!commonGenes.containsKey(m.gene)) {
                commonGenes.put(m.gene,new ArrayList<Mutation>());
                mutCount.put(m.gene, new ArrayList<Mutation>());
            }
            commonGenes.get(m.gene).add(m);

            // For purposes of recessive filtering below
            boolean inTree = false;
            for (Mutation mCount : mutCount.get(m.gene)) {
                if (mCount.equalsOrUnsequenced(m) || mCount.equalsMutation(m)) {
                    inTree = true;
                    break;
                }
            }
            if (!inTree) {
                mutCount.get(m.gene).add(m);
            }
        }

        // if recessive, non-consanguineous, then cut out all non-X genes that have < 2 mutations
        if (flags.RECESSIVE && (!consanguineous)) {

            ArrayList<String> keysToKeep = new ArrayList<String>();
            for (String gene : commonGenes.keySet()) {
                if (commonGenes.get(gene).get(0).chromosome.equalsIgnoreCase("x")) {
                    keysToKeep.add(gene);
                }
                else if (mutCount.get(gene).size() > 1) {
                    // check to make sure each mutation is distinct (not just different transcript)
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
