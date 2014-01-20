package com.example.mutationfilter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.TreeMap;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/11/13
 * Time: 11:56 AM
 * To change this template use File | Settings | File Templates.
 *
 * Data Structure for Storing Centimorgan coordinates contained in genetic_map_GRCh37_chr*.txt files.
 * To be used for homozygous region filtering only.
 *
 * Files are hardcoded in.  Change location of files if running java program elsewhere.
 *
 * Logic: find position closest to candidate position for endpoints of Homozygous region
 * and take difference between cM coordinates for each position
 *
 * Note: Homozygous Region Filtering will be an approximation.
 * Note: chromosome strings stored just as number of character (no "chr")
 *
 */
public class CmToMbMap {

    // Constants
    private static String FILE_DIR = "/Users/kaitlinstouffer/Dropbox/IdeaProjects/MutationFilter/genetic_map_HapMapII_GRCh37/";

    // Instance Variables
    private TreeMap<String,TreeMap<Integer,Double>> rates;

    public CmToMbMap() throws IOException {
        String autosomeFiles = FILE_DIR + "genetic_map_GRCh37_chr";
        rates = new TreeMap<String, TreeMap<Integer,Double>>();
        for (int i = 1; i < 23; i++) {
            // Make new Tree Map for Chromosome
            TreeMap<Integer,Double> tree = new TreeMap<Integer,Double>();

            // Read in file (first line is header and second column is position, fourth column is cM coordinate)
            String f = autosomeFiles + String.valueOf(i) + ".txt";
            TabFileReader reader = new TabFileReader(f);

            String[] line = reader.readColumnsArray();
            line = reader.readColumnsArray(); // skip header row

            while (line != null) {
                tree.put(Integer.parseInt(line[1]), Double.parseDouble(line[3]));
                line = reader.readColumnsArray();
            }
            reader.close();

            rates.put(String.valueOf(i),tree);
        }

        // Add X Chromosome
        /* NOTE: X chromosome is divided into 3 sections, each of which dips down to a 0 cM coordinate.
         * If the difference between cM coordinates is -, then assume is not a continuous region.
         *
         * All of X is stored as X in TreeMap
         */
        String[] xFiles = {autosomeFiles + "X_par1.txt", autosomeFiles + "X.txt", autosomeFiles + "X_par2.txt"};
        TreeMap<Integer,Double> xTree = new TreeMap<Integer,Double>();
        for (int i = 0; i < 3; i++) {
            TabFileReader reader = new TabFileReader(xFiles[i]);
            String[] line = reader.readColumnsArray();
            line = reader.readColumnsArray();

            while (line != null) {
                xTree.put(Integer.parseInt(line[1]), Double.parseDouble(line[3]));
                line = reader.readColumnsArray();
            }
            reader.close();
        }

        rates.put("X", xTree);
    }

    // Returns the cM coordinate for closest position to given location (chromosome and bp position)
    // Chromosome is given as String object to accommodate X chromosome
    public double getcM(String chrom, int bp) {
        TreeMap<Integer, Double> tree = rates.get(chrom);
        // look through set of keys and find one that is closest (LOOK UP HOW TO DO THIS; MIGHT HAVE TO ENSURE KEYS ARE SORTED)
        Integer above = tree.ceilingKey(bp);
        Integer below = tree.floorKey(bp);

        if (above == null) {
            if (below != null) {
                return tree.get(below);
            }
            else {
                System.err.println("ERROR: Base Pair Request in cM to Mb Mapping Not Valid.");
                return -1;
            }
        }
        else if (below == null) {
            return tree.get(above);
        }
        else {
            if (Math.abs(above-bp) <= Math.abs(bp-below)) {
                return tree.get(above);
            }
            else
                return tree.get(below);
        }

    }

}
