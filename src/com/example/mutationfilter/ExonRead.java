package com.example.mutationfilter;

import java.io.FileReader;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.TreeMap;
import javax.script.ScriptContext;
import javax.script.SimpleScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

/**
 * Created by kaitlinstouffer on 1/8/14.
 *
 * Class that calls python script reads_only_locs.py to check whether individual
 * locations in files are sequenced. Alters possibleLocs of FamilyDataGroup if
 * certain locations are not sequenced in those individuals.
 *
 * UPDATE: 8/1/2014 reads_only_locs.py is incorporated into method rather than calling separate script
 */
public class ExonRead {

    // constants
    //public static final String PYTHON_SCRIPT_LOC = "~/Dropbox/IdeaProjects/ExonReads/reads_only_locs.py"; // location of reads_only_loc.py script

    // instance variables
    private TreeMap<Integer,ArrayList<String>> locsToRead;
    private FamilyDataGroup family;
    private String bamFileDir;

    public ExonRead(TreeMap<Integer,ArrayList<String>> ltr, FamilyDataGroup fdg, Flags flags) {
        locsToRead = ltr;
        family = fdg; // reference to family in FilterProgram so can alter possibleLocs directly
        bamFileDir = flags.BAM_DIR;
    }

    /*
     * Iterates over all individuals in a family and checks whether each location
     * in locsToRead for that individual is sequenced.  If not, adds that location
     * as a new "mutation" to the possibleLocs.
     */
    public void checkSequencing() throws ScriptException {

        ScriptEngine engine = new ScriptEngineManager().getEngineByName("python");
        // Start Python Script
        // imports: ensure pysam library is available
        engine.eval("import sys");
        engine.eval("import getopt");
        engine.eval("import pysam");


        // Assume bam file names are the equivalent of .annot.tab files
        for (Integer i : locsToRead.keySet()) {
            // get file name
            int endInd = family.individuals[i].indexOf("annot");
            String bamFileLoc = bamFileDir + family.individuals[i].substring(0,endInd) + "bam"; // assume bamFileDir given with slash

            // establish bamfile and samfile objects
            // TODO: (8/1/2014) ADD EXCEPTIONS IF BAM FILES CAN'T BE OPENED
            engine.put("bamfile", bamFileLoc);
            engine.put("samfile", "");
            engine.eval("samfile = pysam.Samfile( bamfile, \"rb\" )");

            // check format of reference (chr# or #); default is chr#
            engine.put("references", "");
            engine.eval("references = samfile.references");
            engine.put("chromosomeRef", "True");
            engine.eval("chromosomeRef = references[0].startswith('c')");
            // store chromosomeRef in java environment
            boolean chromosomeRef = Boolean.parseBoolean((String) engine.get("chromosomeRef"));

            // iterate over locations to read
            for (String s : locsToRead.get(i)) {
                String[] wholeString = s.split(":");
                engine.put("position", Integer.parseInt(wholeString[1]));
                //engine.put("chromosome", wholeString[0]);
                String chromosome = wholeString[0];

                // Work around mitochondrial and various contig formatting
                if (chromosome.equals("chrM") && !chromosomeRef) {
                    chromosome = "chrMT";
                }

                int indOfGl = wholeString[0].indexOf("gl");
                if (indOfGl >= 0 && !chromosomeRef) {
                    indOfGl += 2;
                    int indOfRand = wholeString[0].indexOf("random");
                    if (indOfRand >= 0) {
                        indOfRand--;
                        chromosome = "chrGL" + chromosome.substring(indOfGl,indOfRand) + ".1";
                    }
                    else {
                        chromosome = "chrGL" + chromosome.substring(indOfGl) + ".1";
                    }

                }
                if (chromosome.split("_").length >= 3) {
                    continue;
                }

                if (!chromosomeRef) {
                    chromosome = chromosome.substring(3);
                }

                String region_string = chromosome + ":" + wholeString[1]
                        + ":" + (Integer.parseInt(wholeString[1]) + 1);
                engine.put("region_string", region_string);
                engine.eval("iter = samfile.pileup(region=region_string,max_depth=20000)");
                engine.put("count", "0");

                engine.eval("count = len([1 for x in iter if (x.pos == (position-1)])");

                // if count = 0, location was not sequenced
                if (Integer.parseInt((String) engine.get("count")) == 0) {
                    family.possibleLocs.get(i).add(new Mutation(s));
                }

            }

        }

    }

}
