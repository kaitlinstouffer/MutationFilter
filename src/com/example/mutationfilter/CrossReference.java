package com.example.mutationfilter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 2:09 PM
 * To change this template use File | Settings | File Templates.
 *
 * Class that compares mutations found in affected individuals to mutations in individuals without these phenotypes.
 * Facilitates the elimination of sequencing errors and artifacts in a list of possible pathogenic mutations.
 *
 * Note: CrossReference Comparison is only made off of location, not comparing whole mutation (5/12/2013)
 * i.e. if two people have a different mutation at same location, then that location is eliminated
 */

public class CrossReference {

    // Constants
    public String EXTENSION = ".annot.tab"; // extension of filenames to consider in referenceDirectories
    // TODO: might want to change this to be an instance variable and flag set

    // Instance Variables
    public String[] referenceDirectories; // list of directories with exomes to compare to affected individuals
    public int refCutOff; // input cutoff for how many individuals must possess mutation to eliminate
    public String[] filesToIgnore;
    public ArrayList<Mutation> possibleLocs;
    public ArrayList<String> possibleLocsString; // Use to compare just locations (not whole mutations)
    public ArrayList<Mutation> filteredLocs; // stores locations kept after cross reference

    // Constructors
    public CrossReference(String[] refDir, int cut, String[] filesIg, ArrayList<Mutation> possLocs) {
        referenceDirectories = refDir;
        refCutOff = cut;
        filesToIgnore = filesIg;
        possibleLocs = possLocs;
        possibleLocsString = new ArrayList<String>();
        for (Mutation m : possibleLocs) {
            // Make list of possible locations (without duplicates)
            if (!possibleLocsString.contains(m.chromLoc)) {
                possibleLocsString.add(m.chromLoc);
            }

        }
    }

    // Methods

    /* Filters for locations that are not contained in any of the files in the
    reference directories.

    @cutoff = specifies how many files in reference directories must contain a
    mutation in order for it to be eliminated from possible locations (default is 1)
     */

    public void filterPossibleLocs() throws IOException {
        int[] removeLoc = new int[possibleLocsString.size()];

        // make list of all of locs to enumerate
        ArrayList<String> filesToEnumerate = new ArrayList<String>();
        // Only keep files that end in .annot.tab
        for (String s : referenceDirectories) {
            File[] files = new File(s).listFiles();
            for (File f : files) {
                if (f.isFile() && f.getName().endsWith(EXTENSION)) {
                    filesToEnumerate.add(s + f.getName());
                }
            }
        }
        // TODO: remove
        System.out.println("in filterPossibleLocs and has made list of files.");
        for (String s : filesToEnumerate) {
            if (Arrays.asList(filesToIgnore).contains(s)) {
                continue;
            }
            System.out.println("reading file " + s);
            TabFileReader reader = new TabFileReader(s);
            String[] cols = reader.readColumnsArray();
            while (cols != null) {
                if (cols.length < 2) {
                    continue;
                }
                if (cols[0].startsWith("c") || cols[0].startsWith("C")) {
                    cols[0] = cols[0].substring(3);
                }
                int loc = possibleLocsString.indexOf(cols[0]);
                if (loc >= 0) {
                    removeLoc[loc]++;
                }
                cols = reader.readColumnsArray();
            }
            reader.close();
        }

        // Remove found locations from list
        ArrayList<String> filteredLocsString = new ArrayList<String>();
        for (int i = 0; i < removeLoc.length; i++) {
            if (removeLoc[i] < refCutOff) {
                filteredLocsString.add(possibleLocsString.get(i));
            }
        }

        // Turn String locations back into mutation objects
        filteredLocs = new ArrayList<Mutation>();
        for (Mutation m : possibleLocs) {
            if (filteredLocsString.contains(m.chromLoc)) {
                filteredLocs.add(m);
            }
        }
        // TODO: remove
        System.out.println("finished filtering locs for cross reference");
    }


}
