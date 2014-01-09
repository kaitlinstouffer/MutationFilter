package com.example.mutationfilter;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 11/25/13
 * Time: 6:31 PM
 * To change this template use File | Settings | File Templates.
 *
 * Data Structure for describing SNPs and Indels found in an exome.
 * TO DO: implements Comparable for an ordering for intersection
 */
public class Mutation {

    // Instance Variables
    public String chromosome; // store chromosome number in either chr# or just #
    public int locationNum; // Location on Chromosome
    public String refAllele; // Reference Base at this Position
    public String altAllele; // Alternative (observed) base at this position
    public String transcript;
    public String chromLoc; // concatenated chromosome and location number
    public double snpFreq; // global frequency of snp from 1000 genomes
    public boolean isSnp; // true if location is snp registered in dbsnp
    public boolean isSequenced; // true if location was sequenced (used only in exonRead comparison)
    public String gene; // used for making gene banks in comparison

    // Constructors
    public Mutation(String c, int l, String ra, String aa, String t, String cl,
                    double sf, boolean s, String g) {
        chromosome = c;
        locationNum = l;
        refAllele = ra;
        altAllele = aa;
        transcript = t;
        chromLoc = cl;
        snpFreq = sf;
        isSnp = s;
        isSequenced = true;
        gene = g;

    }

    // Make Mutation object for location that is not sequenced in a file
    // Always stores chromosome in format of number rather than chr#
    public Mutation(String c) {
        if (c.startsWith("c")) {
            chromLoc = c.substring(3);
        }
        String[] loc = chromLoc.split(":");
        chromosome = loc[0];
        locationNum = Integer.parseInt(loc[1]);

        transcript = null;
        refAllele = null;
        altAllele = null;
        snpFreq = -1;
        isSnp = false;
        isSequenced = false;
        gene = null;
    }

    // Set Mutation object based on header indices and row of tab file
    public Mutation(Flags f, String[] row) {
        String[] loc = row[0].split(":");
        chromLoc = row[0];
        chromosome = loc[0];
        locationNum = Integer.parseInt(loc[1]);

        transcript = row[f.TRANSCRIPT];
        String[] change = row[f.REF].split(">");
        refAllele = change[0];
        altAllele = change[1];
        if (row[f.SNP].contains("rs")) {
            isSnp = true;
        }
        else if (!f.RS_ONLY && !row[f.SNP].equals(".")) {
            isSnp = true;
        }
        else {
            isSnp = false;
            snpFreq = -1;
        }

        //TO DO (12/12)....FINISH THIS METHOD
        if (f.GMAF > 0 && isSnp) {
            if (row[f.GMAF] != null && !row[f.GMAF].isEmpty()) {
                String[] ss = row[f.GMAF].split("&");
                boolean freqSet = false;
                for (String s : ss) {
                    if (s.startsWith(altAllele)) {
                        snpFreq = Double.parseDouble(s.split(":")[1]);
                        freqSet = true;
                        break;
                    }
                }
                // subtract all of frequencies from 1 if not set
                if (!freqSet) {
                    snpFreq = 1;
                    for (String s : ss) {
                        snpFreq -= Double.parseDouble(s.split(":")[1]);
                    }

                }
            }
            else {
                snpFreq = 0; // set frequency to 0 if not given
            }
        }

        gene = row[f.GENE];
    }

    // Equals Methods
    public boolean equals(Mutation m) {
        // check if are same object
        if (m == this)
            return true;
        if (this.chromosome.equals(m.chromosome) && (this.locationNum == m.locationNum)
                && this.refAllele.equals(m.refAllele) && this.altAllele.equals(m.altAllele)
                && this.transcript.equals(m.transcript))
            return true;
        return false;
    }

    // Returns true if Mutation objects describe same base change at same location
    public boolean equalsMutation(Mutation m) {
        if (m == this)
            return true;
        if (this.chromosome.equals(m.chromosome) && (this.locationNum == m.locationNum)
                && this.refAllele.equals(m.refAllele) && this.altAllele.equals(m.altAllele))
            return true;
        return false;

    }

    // Returns true if Mutation objects have same location but one is not sequenced
    // or if regular equals method is true
    public boolean equalsOrUnsequenced(Mutation m) {
        if (m.isSequenced && this.isSequenced) {
            return this.equals(m);
        }
        // at least one of objects is not sequenced, so chromosomes and locs only must be the same
        else {
            if (this.chromosome.equals(m.chromosome) && this.locationNum == m.locationNum) {
                return true;
            }
            else {
                return false;
            }
        }
    }

}
