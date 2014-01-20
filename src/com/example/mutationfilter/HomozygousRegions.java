package com.example.mutationfilter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/10/13
 * Time: 2:31 PM
 * To change this template use File | Settings | File Templates.
 *
 * Class that documents and stores coordinates of all regions in exome that are a given percentage
 * and given length (in cM) homozygous.  Assumes Homozygous Regions always start and end with Homo.
 *
 * Built only based off of locations in exomes.
 *
 * For use with Consanguineous Homozygous Filtering.
 */
public class HomozygousRegions {

    // Instance Variables
    private double percentage = 0.0;
    private int cMLength = 0;  // should be formed from commandline arg
    private TreeMap<String,ArrayList<EndPoints>> homoRegions; // chromosome Number and arraylist of pairs of start and end coordinates
    private Flags flags;

    private class EndPoints {
        // instance variables
        private int start;
        private int end;

        public EndPoints(int s, int e) {
            start = s;
            end = e;
        }

        // Returns true if i is contained in the interval between s and e (includes endpoints)
        public boolean contains(int i) {
            if (i <= end && i >= start)
                return true;
            else
                return false;

        }

    }

    // Returns true if the ArrayList of EndPoints corresponding to given chromosome
    // contains an interval that the given integer is contained in
    // don't assume ArrayList is sorted necessarily
    public boolean containsLoc(String chrom, int i) {
        ArrayList<EndPoints> a = homoRegions.get(chrom);
        for (EndPoints e : a) {
            if (e.contains(i)) {
                return true;
            }
        }
        return false;
    }

    // Constructor of Homozygous Regions object (made from individual file)
    public HomozygousRegions(double per, int len, Flags f, String filename, CmToMbMap conversion) throws IOException {
        percentage = per;
        cMLength = len;
        flags = f;
        homoRegions = new TreeMap<String,ArrayList<EndPoints>>();

        // if HOMO header not found, cannot make Homozygous Regions
        if (flags.HOMO < 0) {
            System.err.println("ERROR: HOMOZYGOSITY column not found in files or not set.");
            return;
        }

        // read through file (assume it is sorted)
        TabFileReader reader = new TabFileReader(filename);
        String[] line = reader.readColumnsArray();
        // skip header and extra information lines
        while (!line[0].contains(":") || (line.length < 4)) {
            line = reader.readColumnsArray();
        }

        // Approach 2:
        /* Iterate over lines in file and keep track of blocks of at least cMLength.
         * If a block of cMLength is 95% homozygous, then put into list.  Otherwise,
         * slide down one location and re-calculate.
         */
        while (line != null) {
            // start finding cM block
            ArrayList<String[]> linesOfBlock = new ArrayList<String[]>();
            String chromLoc = line[0]; // don't count the same location multiple times
            String[] locParts = chromLoc.split(":");
            int startOfRegion = Integer.parseInt(locParts[1]);
            String chromosome = locParts[0];
            double startOfRegionCM = conversion.getcM(chromosome, startOfRegion);
            int endOfRegion = startOfRegion;
            double endOfRegionCM = startOfRegionCM;
            int numOfHet = 0;

            // walk through locations to find until reach cMLength
            while ((endOfRegion >= startOfRegion) && ((endOfRegionCM - startOfRegionCM) < cMLength)) {
                if (line[flags.HOMO].contains("HET")) {
                    numOfHet++;
                }
                linesOfBlock.add(line);
                while (line != null && line[0].equals(chromLoc)) {
                    line = reader.readColumnsArray();
                }
                if (line == null)
                    break;
                chromLoc = line[0];
                endOfRegion = Integer.parseInt(line[0].split(":")[1]);
                endOfRegionCM = conversion.getcM(chromosome, endOfRegion);

            }
            boolean changeChrom = false;
            // if chromosomes changed, while finding this last one, then skip adding this block
            if (endOfRegion <= startOfRegion || line == null) {
                // check last element that was added
                String[] lastLoc = linesOfBlock.get(linesOfBlock.size()-1)[0].split(":");
                if ((conversion.getcM(lastLoc[0], Integer.parseInt(lastLoc[1])) - startOfRegionCM) >= cMLength) {
                    // keep this block if percentage is correct
                    if ((1.0 - (numOfHet / (double) linesOfBlock.size())) >= percentage) {
                        if (!homoRegions.containsKey(chromosome)) {
                            homoRegions.put(chromosome, new ArrayList<EndPoints>());
                        }
                        homoRegions.get(chromosome).add(new EndPoints(startOfRegion,Integer.parseInt(lastLoc[1])));
                        // TODO: remove print debugging statement
                        System.out.println("added " + chromosome + ":" + startOfRegion + "-" + lastLoc[1]);
                    }

                }

            }


            else {
                // slide region down chromosome until is 95% homozygous
                if (line[flags.HOMO].contains("HET")) {
                    numOfHet++;
                }
                linesOfBlock.add(line);
                while ((1.0 - (numOfHet / (double) linesOfBlock.size())) < percentage) {
                    // remove first element and reset starting coordinates
                    linesOfBlock.remove(0);
                    String[] newStart = linesOfBlock.get(0)[0].split(":");
                    startOfRegion = Integer.parseInt(newStart[1]);
                    startOfRegionCM = conversion.getcM(chromosome, startOfRegion);

                    // add coordinates until reach large enough length
                    while ((endOfRegionCM - startOfRegionCM) < cMLength) {
                        line = reader.readColumnsArray();
                        if (line == null) {
                            break;
                        }
                        if (line.equals(chromLoc)) {
                            continue;
                        }
                        String[] loc = line[0].split(":");
                        if (!loc[0].equals(chromosome)) {
                            changeChrom = true;
                            break;
                        }
                        chromLoc = line[0];
                        endOfRegion = Integer.parseInt(loc[1]);
                        endOfRegionCM = conversion.getcM(chromosome, endOfRegion);
                        if (line[flags.HOMO].contains("HET")) {
                            numOfHet++;
                        }
                        linesOfBlock.add(line);

                    }
                    if (line == null || changeChrom) {
                        break;
                    }
                }
                if (line == null || changeChrom) {
                    continue;
                }
                else {
                    // add block to list and get new line
                    if (!homoRegions.containsKey(chromosome)) {
                        homoRegions.put(chromosome,new ArrayList<EndPoints>());
                    }
                    homoRegions.get(chromosome).add(new EndPoints(startOfRegion, endOfRegion));
                    System.out.println("added " + chromosome + ":" + startOfRegion + "-" + endOfRegion);
                    line = reader.readColumnsArray();
                    while (line != null && line[0].equals(chromLoc)) {
                        line = reader.readColumnsArray();
                    }
                }
            }

        }

        // Approach 1
        /*
        // iterate over lines in file and keep track of start and end of homozygous regions
        int startOfRegion = 0;
        int endOfRegion = 0;
        int numOfLocs = 0;
        int numOfHet = 0;
        String chromosome = null;
        boolean inHomo = false;
        while (line != null) {
            // if not in Homozygous region already, look for first homozygous
            if (!inHomo) {
                if (line[flags.HOMO].contains("HOMO")) {
                    inHomo = true;
                    String[] locParts = line[0].split(":");
                    startOfRegion = Integer.parseInt(locParts[1]);
                    chromosome = locParts[0];
                    endOfRegion = Integer.parseInt(locParts[1]);
                    numOfLocs = 1;
                }
            }
            else {
                // if Homo, then update endOfRegion and numOfLocs
                // check for change in chromosome or region coordinates being overturned
                if (line[flags.HOMO].contains("HOMO")) {
                    String[] locParts = line[0].split(":");
                    if (chromosome.equals(locParts[0])) {
                        int newEndOfRegion = Integer.parseInt(locParts[1]);
                        if (conversion.getcM(chromosome,newEndOfRegion) - conversion.getcM(chromosome,startOfRegion) > 0) {
                            endOfRegion = newEndOfRegion;
                            numOfLocs++;
                        }
                        // if cM region now reads negative, then is a jump in the coordinates, so end homozygous region
                        else {
                            checkHomoRegion(startOfRegion,endOfRegion,chromosome,numOfLocs,numOfHet,conversion);
                            numOfLocs = 1;
                            numOfHet = 0;
                            startOfRegion = newEndOfRegion;
                            endOfRegion = newEndOfRegion;
                        }

                    }
                    else {
                        checkHomoRegion(startOfRegion,endOfRegion,chromosome,numOfLocs,numOfHet,conversion);
                        numOfLocs=1;
                        numOfHet = 0;
                        startOfRegion = Integer.parseInt(locParts[1]);
                        endOfRegion = startOfRegion;
                    }
                }
                else {
                    String[] locParts = line[0].split(":");
                    int newEndOfRegion = Integer.parseInt(locParts[1]);

                    // if jump chromosomes, check for putting last region in TreeMap and reset
                    if (!locParts[0].equals(chromosome)) {
                        checkHomoRegion(startOfRegion,endOfRegion,chromosome,numOfLocs,numOfHet,conversion);
                        inHomo = false;
                        numOfLocs = 0;
                        numOfHet = 0;
                    }
                    // if percentage is still acceptable, then keep retaining same Homozygous region
                    if ((double)(numOfLocs+1 - (numOfHet + 1))/(numOfLocs+1) >= percentage) {
                        endOfRegion = newEndOfRegion;
                        numOfHet++;
                        numOfLocs++;
                    }

                    else {
                        checkHomoRegion(startOfRegion,endOfRegion,chromosome,numOfLocs,numOfHet,conversion);
                        inHomo = false;
                        numOfLocs = 0;
                        numOfHet = 0;
                    }
                }

            }
            line = reader.readColumnsArray();
        }
        */
        reader.close();


    }

    // Helper Method to check if region should be added as homozygous region
    private void checkHomoRegion(int start, int end, String chrom, int numLocs, int numHet, CmToMbMap con) {
        double startCM = con.getcM(chrom, start);
        double endCM = con.getcM(chrom, end);
        if (((1.0 - (numHet/(double)numLocs)) >= percentage) && ((endCM - startCM) >= cMLength)) {
            // make Homozygous Region
            EndPoints ep = new EndPoints(start,end);
            if (homoRegions.containsKey(chrom)) {
                homoRegions.get(chrom).add(ep);
            }
            else {
                homoRegions.put(chrom, new ArrayList<EndPoints>());
                homoRegions.get(chrom).add(ep);
            }
        }

    }
}
