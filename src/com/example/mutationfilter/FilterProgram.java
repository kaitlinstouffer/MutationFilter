package com.example.mutationfilter;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import static com.example.mutationfilter.FamilyDataGroup.makeDataGroups;
import static com.example.mutationfilter.FamilyDataGroup.setBamFileNames;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 2:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class FilterProgram {

    // Outline for Program (6/12/2013)
    /*
     * Pass the tab-separated row through each filter
     * If passes all of filters, then add as a mutation to list
     * No need to do intersection per individual, just have final list of all of mutations/locations that passed all of filters
     *
     * Then before compare individuals, look at all locations and see if were read or not using reads_only_locs.py
     * on telemachus; if not read in a person's file, then add to their list of possible locations, but with indication
     * that were not read
     *
     * Do gene banks / organization adn comparing across families
     *
     * Do cross reference
     *
     *
     */

    public static void main(String[] args) throws IOException, ScriptException {
        // Create Family Data Group
        ArrayList<FamilyDataGroup> families = makeDataGroups(args[0]);

        // args[1] denotes root for where to write output and temporary files
        String outputRoot = args[1];

        // Parse Commandline Flags to make list of Filters to Run
        Flags flags = new Flags();
        flags.parseCommandLineArgs(args);

        // set Consanguinity flags per family
        FamilyDataGroup.setConsanguinity(families, flags);

        ArrayList<Filter> filters = new ArrayList<Filter>();
        // Check for settings of each filter and create filter if specified
        if (flags.FILTER_INHERITANCE) {
            InheritanceFilter inF = new InheritanceFilter(flags);
            filters.add(inF);
        }
        if (flags.FILTER_MUT_TYPE) {
            MutationTypeFilter mF = new MutationTypeFilter(flags);
            filters.add(mF);
        }
        if (flags.FILTER_QUALITY) {
            QualityFilter qF = new QualityFilter(flags);
            filters.add(qF);
        }
        if (flags.FILTER_SNP) {
            SnpFilter sF = new SnpFilter(flags);
            filters.add(sF);
        }
        // TODO: REMOVE
        System.out.println("Filters = " + filters.size());

        // Parse through files and find candidate locations for each individual
        // NOTE: Make all comparisons with location as number/character:loc, no "chr"; trim before running through filters
        // but save information of whether file contains chr references or non-chr references in family data group
        for (FamilyDataGroup family : families)
        {
            ArrayList<ArrayList<Mutation>> familyMuts = new ArrayList<ArrayList<Mutation>>();
            family.chrReferences = new boolean[family.individuals.length];
            int count = 0;
            for (String ind : family.individuals) {
                ArrayList<Mutation> indMuts = new ArrayList<Mutation>();
                // if family is consanguineous and filtering based on homozygous regions, make new homozygous region for individual
                // Inheritance filter will be first in list of filters
                if (flags.FILTER_INHERITANCE && ((flags.CONSANGUINEOUS > 0) && family.consanguineous)) {
                    InheritanceFilter inF = (InheritanceFilter) filters.get(0);
                    inF.setHomoRegions(ind);
                }

                TabFileReader reader = new TabFileReader(ind);
                String[] line = reader.readColumnsArray();
                // skip header row
                while (!line[0].contains(":") || (line.length < 4)) {
                    if (line[0].contains("osition")) {
                        flags.parseHeader(line);
                        // update all of filters
                        for (Filter f : filters) {
                            f.updateFlagHeaders(line);
                        }
                    }
                    line = reader.readColumnsArray();
                }

                // Parse actual locations
                while (line != null) {
                    // remove "chr" from location
                    if (line[0].startsWith("c")) {
                        line[0] = line[0].substring(3);
                        family.chrReferences[count] = true;
                    }
                    boolean passFilters = true;
                    for (Filter f : filters) {
                        if (!f.pass(line, family)) {
                            passFilters = false;
                            break;
                        }
                    }

                    // if passed all of filters, make location and add to list
                    if (passFilters) {
                        Mutation m = new Mutation(flags, line);
                        indMuts.add(m);
                    }
                    line = reader.readColumnsArray();
                }
                reader.close();
                // Add individual mutation list to family mutation list
                familyMuts.add(indMuts);
                count++;
            }
            // set familyMuts as family data group instance variable
            family.possibleLocs = familyMuts;
            // TODO: REMOVE
            for (ArrayList<Mutation> am : family.possibleLocs) {
                System.out.println("number of possible locs is " + am.size());
            }
        }

        // If Exon Reads Used, make list of possible locations for each family to check if read
        if (flags.NUM_EXON_READ > 0) {
            // set bam files in families
            setBamFileNames(families,flags.BAM_FILE_NAMES);
            for (FamilyDataGroup fdg : families) {
                // Only do Exon Reads if number of individuals that must have mutation in family is not = everyone
                if (flags.NUM_EXON_READ < fdg.individuals.length) {
                    // Integer key represents number of individual in family data group
                    TreeMap<Integer,ArrayList<String>> locsToRead = fdg.getLocsToRead(flags.NUM_EXON_READ);
                    // check if locs are read for each individual and if so, add to possibleLocs
                    if (locsToRead != null) {
                        ExonRead er = new ExonRead(locsToRead, fdg);
                        er.checkSequencing();
                    }
                }
            }
        }

        // Make list of possible mutations assuming each individual in family has mutation
        // in his or her list of possibleLocs

        // Make list of genes associated with possible mutations
        // Recessive, non-consanguineous families must have 2 mutation in same gene to be considered
        for (FamilyDataGroup fdg : families) {
            fdg.setCommonLocs();
            fdg.setCommonGenes(flags);
        }

        // Print Potentially Pathogenic Genes in Families: all Printing will be grouped by gene
        // over family unless specified in output flag
        String outputFile = outputRoot + "_Genes.txt";
        if (flags.PRINT_BY_FAM) {

            printMutationsByGeneFamSeparate(flags, families, outputFile);
        }
        else {
            System.out.println("in print");
            printMutationsByGene(flags, families, outputFile);
            System.out.println("out of print");
        }

        if (flags.CROSS_REF_DIR != null) {
            // Complete Cross Reference to Eliminate Artifacts in Sequencing
            for (FamilyDataGroup fdg : families) {
                CrossReference cr = new CrossReference(flags.CROSS_REF_DIR,flags.REF_CUTOFF,fdg.individuals,fdg.commonLocs);
                cr.filterPossibleLocs();
                fdg.commonLocs = cr.filteredLocs;
                // TODO: REMOVE PRINT STATEMENTS
                System.out.println("Filtered locs in cross reference.");
                fdg.setCommonGenes(flags); // reset gene dictionary with new subset of locs
                System.out.println("Reset locs for family.");
            }

            outputFile = outputRoot + "_CrossRef_Genes.txt";
            if (flags.PRINT_BY_FAM) {
                printMutationsByGeneFamSeparate(flags, families, outputFile);
            }
            else {
                printMutationsByGene(flags,families,outputFile);
            }
        }

        // TODO: implement ability to filter out mutations or rank by how many families have mutations in those genes

    }

    // Print Methods

    /* Prints the candidate mutations for each family, grouped by gene.
     * Reads in exact row in input file and replicates it to ensure all data elements
     * are reiterated.
     *
     * Families are printed separately.
     */
    public static void printMutationsByGeneFamSeparate(Flags flags, ArrayList<FamilyDataGroup> families, String outputFile) throws IOException {
        TabFileWriter writer = new TabFileWriter(outputFile);
        int count = 0;
        for (FamilyDataGroup fdg : families) {
            writer.write("Family " + count + ":\n");

            /*
            for (String indName : fdg.individuals) {
                writer.write(indName + " ");
            }
            */
            writer.write("\n");
            writer.write("--------------------------------------------------------------------------------------\n");
            writer.write("--------------------------------------------------------------------------------------\n");
            TreeMap<String, TreeMap<Mutation,ArrayList<ArrayList<String[]>>>> linesToPrint = getInputGeneLines(fdg, flags);
            for (String gene : linesToPrint.keySet()) {
                for (Mutation m : linesToPrint.get(gene).keySet()) {
                    int indvCount = 0;
                    writer.write("* " + m.chromLoc + " *\n");
                    for (ArrayList<String[]> as : linesToPrint.get(gene).get(m)) {
                        writer.write(fdg.individuals[indvCount] + "\n");
                        for (String[] s : as) {
                            writer.writeColumns(s);
                        }
                        writer.write("\n");
                        indvCount++;
                    }
                }
                writer.write("--------------------------------------------------------------------------------------\n");
            }
            writer.write("--------------------------------------------------------------------------------------\n");
            writer.write("--------------------------------------------------------------------------------------\n");
            writer.write("\n");
            writer.write("\n");
            count++;
        }
        writer.close();
    }

    /* Prints the candidate mutations for each family, grouped by gene
     * over grouping by family.
     */
    public static void printMutationsByGene(Flags flags, ArrayList<FamilyDataGroup> families, String outputFile) throws IOException {
        TabFileWriter writer = new TabFileWriter(outputFile);
        ArrayList<TreeMap<String, TreeMap<Mutation,ArrayList<ArrayList<String[]>>>>> linesToPrint =
                new ArrayList<TreeMap<String, TreeMap<Mutation,ArrayList<ArrayList<String[]>>>>>();
        ArrayList<String> allGenes = new ArrayList<String>();
        for (FamilyDataGroup fdg : families) {
            linesToPrint.add(getInputGeneLines(fdg,flags));
            for (String gene : fdg.commonGenes.keySet()) {
                if (!allGenes.contains(gene))
                    allGenes.add(gene);
            }
        }

        // Iterate over all genes and print out family data if exists for given gene
        for (String gene : allGenes) {
            if (flags.PRINT_SHORT_GENE) {
                int total = 0;
                for (int k = 0; k < families.size(); k++) {
                    TreeMap<Mutation,ArrayList<ArrayList<String[]>>> temp = linesToPrint.get(k).get(gene);
                    if (temp != null) {
                        total++;
                    }
                }
                if (total < 2) {
                    continue;
                }
            }
            writer.write(gene + ":\n");
            writer.write("--------------------------------------------------------------------------------------\n");
            writer.write("--------------------------------------------------------------------------------------\n");

            for (int i = 0; i < families.size(); i++) {
                TreeMap<Mutation,ArrayList<ArrayList<String[]>>> famLines = linesToPrint.get(i).get(gene);
                if (famLines == null) {
                    continue;
                }

                writer.write("Family " + i + ":\n");
                /*
                for (String indName : families.get(i).individuals) {
                    writer.write(indName + " ");
                }
                */
                writer.write("\n");
                for (Mutation m : famLines.keySet()) {
                    writer.write("* " + m.chromLoc + " *\n");
                    int indvCount = 0;
                    for (ArrayList<String[]> as : famLines.get(m)) {
                        writer.write(families.get(i).individuals[indvCount] + "\n");
                        for (String[] s : as) {
                            writer.writeColumns(s);
                        }
                        writer.write("\n");
                        indvCount++;
                    }
                }
                writer.write("--------------------------------------------------------------------------------------\n");
                writer.write("\n");
            }
            writer.write("--------------------------------------------------------------------------------------\n");
            writer.write("--------------------------------------------------------------------------------------\n");
            writer.write("\n");
        }
        writer.close();

    }

    /*
     * Reads in input files, one at a time, and registers the rows that correspond to the
     * locations given as possible in common genes.
     *
     * Treats mutations with same change and at same locations as the same.
     */
    private static TreeMap<String, TreeMap<Mutation,ArrayList<ArrayList<String[]>>>> getInputGeneLines(FamilyDataGroup fdg, Flags flags) throws IOException {
        TreeMap<String, TreeMap<Mutation, ArrayList<ArrayList<String[]>>>> inputLines =
                new TreeMap<String, TreeMap<Mutation, ArrayList<ArrayList<String[]>>>>();
        // Make slot for each gene
        for (String gene : fdg.commonGenes.keySet()) {
            inputLines.put(gene,new TreeMap<Mutation, ArrayList<ArrayList<String[]>>>());
            // add Mutation for each gene and ArrayList per individual
            for (Mutation m : fdg.commonGenes.get(gene)) {
                inputLines.get(gene).put(m,new ArrayList<ArrayList<String[]>>());
                for (int i = 0; i < fdg.individuals.length; i++) {
                    inputLines.get(gene).get(m).add(new ArrayList<String[]>());
                }
            }
        }
        for (int i = 0; i < fdg.individuals.length; i++) {
            String file = fdg.individuals[i];
            TabFileReader reader = new TabFileReader(file);
            String[] input = reader.readColumnsArray();
            while (input != null) {
                if (inputLines.containsKey(input[flags.GENE])) {
                    boolean changedInput = false;
                    if (input[0].startsWith("c")) {
                        input[0] = input[0].substring(3);
                        changedInput = true;
                    }
                    Mutation m = new Mutation(flags,input);
                    if (inputLines.get(input[flags.GENE]).keySet().contains(m)) {
                        if (changedInput) {
                            input[0] = "chr" + input[0];
                        }
                        // Avoid reference changing data stored in inputLines
                        String[] temp = input;
                        inputLines.get(input[flags.GENE]).get(m).get(i).add(temp);
                    }

                }
                input = reader.readColumnsArray();
            }
            reader.close();
            // If any of lists of genes do not have input line, assume not read and indicate this
            for (String s : inputLines.keySet()) {
                for (Mutation m : inputLines.get(s).keySet()) {
                    if (inputLines.get(s).get(m).get(i).isEmpty()) {
                        String[] toAdd = ("LOCATION\tNOT\tSEQUENCED\tIN\tEXOME").split("\t");
                        inputLines.get(s).get(m).get(i).add(toAdd);
                    }
                }
            }
        }

        // group together lines per person corresponding to mutations that just differ by transcript
        TreeMap<String, ArrayList<ArrayList<Mutation>>> indsToCombine = new TreeMap<String, ArrayList<ArrayList<Mutation>>>();
        for (String s : inputLines.keySet()) {
            indsToCombine.put(s,new ArrayList<ArrayList<Mutation>>());
            // mutations in sorted order
            Mutation mTemp = null;
            ArrayList<Mutation> ai = new ArrayList<Mutation>();
            for (Mutation m : inputLines.get(s).keySet()) {
                if (!m.equalsMutation(mTemp) && !ai.isEmpty()) {
                    indsToCombine.get(s).add(ai);
                    ai = new ArrayList<Mutation>();
                }
                ai.add(m);
                mTemp = m;
            }
            // add last list
            indsToCombine.get(s).add(ai);
        }

            // Combine ArrayList<String[]> per individual corresponding to same mutation,
            // with different transcript only
            for (String s : inputLines.keySet()) {
                TreeMap<Mutation,ArrayList<ArrayList<String[]>>> newVal = new TreeMap<Mutation, ArrayList<ArrayList<String[]>>>();
                for (ArrayList<Mutation> ai : indsToCombine.get(s)) {
                    // Actual Mutation object just referenced by the first one in list
                    Mutation ref = ai.get(0);
                    newVal.put(ref,inputLines.get(s).get(ref));
                    for (int k = 1; k < ai.size(); k++) {
                        ArrayList<ArrayList<String[]>> temp = inputLines.get(s).get(ai.get(k));
                        int j = 0;
                        for (ArrayList<String[]> t : temp) {
                            for (String[] st : t) {
                                newVal.get(ref).get(j).add(st);
                            }
                            j++;
                        }
                    }
                }
                inputLines.put(s, newVal);
            }

        return inputLines;
    }

    // Helper Methods

    /* Returns an ArrayList of locations only (in form of chr:loc) that are
     * represented by ArrayList of Mutation objects.
     *
     * if chrRef is true for given mut, then loc is written with chr reference.
     */
    public static ArrayList<String> getLocs(ArrayList<Mutation> muts, boolean[] chrRef) {
        ArrayList<String> chrLocs = new ArrayList<String>();
        int count = 0;
        for (Mutation m : muts) {
            if (!chrRef[count]) {
                chrLocs.add(m.chromLoc);
            }
            else {
                String s = "chr" + m.chromLoc;
                chrLocs.add(s);
            }
            count++;
        }
        return chrLocs;
    }
}
