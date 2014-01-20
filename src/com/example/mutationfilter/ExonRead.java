package com.example.mutationfilter;

import java.io.*;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Properties;
import javax.script.ScriptContext;
import javax.script.SimpleScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import javax.script.ScriptEngineFactory;
import org.python.core.Py;
import org.python.core.PySystemState;
import org.python.util.PythonInterpreter;
import org.python.core.PyString;

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
    public static final String PYTHON_SCRIPT_LOC = "/Users/kaitlinstouffer/Dropbox/IdeaProjects/MutationFilter/ExonReads/reads_only_locs2.py"; // location of reads_only_loc.py script

    public static final String PYTHONPATH = "/Users/kaitlinstouffer/Library/Python/2.7/lib/python/site-packages/distribute-0.6.34-py2.7.egg:/Users/kaitlinstouffer/Library/Python/2.7/lib/python/site-packages/pysam-0.7.5-py2.7-macosx-10.9-intel.egg:/Users/kaitlinstouffer/Library/Python/2.7/lib/python/site-packages/Cython-0.19.2-py2.7-macosx-10.9-intel.egg:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages:/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old:/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload:/Users/kaitlinstouffer/Library/Python/2.7/lib/python/site-packages:/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/PyObjC:/Library/Python/2.7/site-packages";
    public static final String LOCATIONFILE = "/Users/kaitlinstouffer/Dropbox/IdeaProjects/MutationFilter/bamLocs.txt";
    // instance variables
    private TreeMap<Integer,ArrayList<String>> locsToRead;
    private FamilyDataGroup family;

    public ExonRead(TreeMap<Integer,ArrayList<String>> ltr, FamilyDataGroup fdg) {
        locsToRead = ltr;
        family = fdg; // reference to family in FilterProgram so can alter possibleLocs directly
    }

    /*
     * Iterates over all individuals in a family and checks whether each location
     * in locsToRead for that individual is sequenced.  If not, adds that location
     * as a new "mutation" to the possibleLocs.
     */
    public void checkSequencing() throws ScriptException, IOException {

        // Assume bam file names are the equivalent of .annot.tab files
        for (Integer i : locsToRead.keySet()) {
            // get file name
            String bamFileLoc = family.bamIndividuals[i];
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(LOCATIONFILE)));
            int count = 0;
            for (String s : locsToRead.get(i)) {
                writer.write(s);
                writer.write("\n");
                count++;
            }
            writer.close();
            System.out.println("Checked " + count + " locs");
            System.out.println("size of poss locs before is " + family.possibleLocs.get(i).size());

            // From: http://programmersheaven.com/discussion/415726/invoking-python-script-from-java
            try {
                    Runtime r = Runtime.getRuntime();
                    Process p = r.exec("python " + PYTHON_SCRIPT_LOC + " " + bamFileLoc + " " + LOCATIONFILE);
                    BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    BufferedReader brError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                    p.waitFor();
                    String line = "";
                    while (brError.ready()) {
                        System.out.println(brError.readLine());
                    }
                // add output locs to
                    while (br.ready()) {
                        String loc = br.readLine();
                        if (!loc.isEmpty()) {
                            family.possibleLocs.get(i).add(new Mutation(loc));
                            System.out.println("Added mutation to " + i);
                        }
                    }
                p.destroy();

                }
                catch (Exception e)
                {
                    String cause = e.getMessage();
                    if (cause.equals("python: not found"))
                        System.out.println("No python interpreter found.");
                }
            System.out.println("size of poss locs after is " + family.possibleLocs.get(i).size());
            }

        /*
        // TODO: other possible addition
        // http://sourceforge.net/mailarchive/forum.php?thread_name=18d1bad30805271911l52d7efe9q22f916421b20dfb%40mail.gmail.com&forum_name=jython-users
        //Properties props = new Properties();
        //props.setProperty("python.path", "/Users/kaitlinstouffer/Dropbox/IdeaProjects/MutationFilter/pysam/");
        //PythonInterpreter.initialize(System.getProperties(), props, new String[] {""});

        // TODO: added to see if works by adding where pysam should be located
        // from: https://wiki.python.org/jython/UserGuide#using-jsr-223
        //PySystemState engineSys = new PySystemState();
        //System.out.println(Py.getSystemState());
        //engineSys.path.append(Py.newString("/Users/kaitlinstouffer/Library/Python/2.7/lib/python/site-packages/pysam-0.7.5-py2.7-macosx-10.9-intel.egg/"));
        //engineSys.path.append(Py.newString("/Users/kaitlinstouffer/Dropbox/IdeaProjects/MutationFilter/pysam/"));
        //Py.setSystemState(engineSys);

        PythonInterpreter interp = new PythonInterpreter(null, new PySystemState());

        PySystemState sys = Py.getSystemState();
        sys.path.append(new PyString("/Users/kaitlinstouffer/jython2.2.1/"));
        sys.path.append(new PyString("/Users/kaitlinstouffer/Library/Python/2.7/lib/python/site-packages/pysam-0.7.5-py2.7-macosx-10.9-intel.egg/"));
        ScriptEngine engine = new ScriptEngineManager().getEngineByName("python");
        */
        // TODO: remove (solution #3)
        /*
        try
        {
            Properties props = new Properties();
            props.setProperty("python.path", PYTHONPATH);

            PythonInterpreter.initialize(System.getProperties(), props, new String[0]);
            PythonInterpreter engine = new PythonInterpreter();
            engine.exec("import sys");
            engine.exec("print sys.path");
            engine.exec("import getopt");
            engine.exec("import pysam");
            System.out.println("completed import");
            engine.set("bamfile", family.bamIndividuals[0]);
            System.out.println("completed bamfile statement");
            System.out.println("bamfile is " + engine.get("bamfile"));

        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        */

        // TODO: remove comment
        /*
        // Start Python Script
        // imports: ensure pysam library is available
        // TODO: change back to engine.eval
        engine.eval("import sys");
        engine.eval("import getopt");
        engine.eval("import pysam");



        // Assume bam file names are the equivalent of .annot.tab files
        for (Integer i : locsToRead.keySet()) {
            // get file name
            //int endInd = family.individuals[i].indexOf("annot");
            //String bamFileLoc = bamFileDir + family.individuals[i].substring(0,endInd-1) + bamFileExt; // assume bamFileDir given with slash
            String bamFileLoc = family.bamIndividuals[i];
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
        */

    }

}
