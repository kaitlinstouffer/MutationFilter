package com.example.mutationfilter;

import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 3:56 PM
 * To change this template use File | Settings | File Templates.
 *
 * Class that stores all of the commandline flags and information about input files.
 * Flags divided into Flags regarding how to run filters and which filters to run.
 *
 * NOTES:
 * 1) commandline parsing may malfunction if the flags that require secondary arguments because might
 *      parse an arg that corresponds to next flag rather than same one
 * 2) TODO: consider adding ignoreX flag to allow individual to completely eliminate any sort of X mutation
 *
 */
public class Flags {

    // Instance Variables

    // Header Values: if has -1 as the index for a header string, then not in file
    //POS, DEPTH, PASS, SNP, GENE, HOMO, MUT_TYPE, TRANSCRIPT, GMAF, REF
    public int POS = -1;
    public int DEPTH = -1;
    public int PASS = -1;
    public int SNP = -1;
    public int GENE = -1;
    public int HOMO = -1;
    public int MUT_TYPE = -1;
    public int TRANSCRIPT = -1;
    public int GMAF = -1;
    public int REF = -1;

    // Inheritance Flags
    // X_LINKED, RECESSIVE, (INCLUDE_HOMO)
    public boolean X_LINKED = false;
    public boolean RECESSIVE = false;
    public boolean INCLUDE_HOMO = false; // might be deprecated

    //Filter Specific Flags
    public int CONSANGUINEOUS = -1;
    public ArrayList<Integer> CONSANGUINEOUS_FAMILIES = null;
    public double SNP_FREQ = 0.0;
    public String[] CROSS_REF_DIR;
    public int DEPTH_CUTOFF = -1; // if none given, then defaults to 10 for non-consanguineous and 20 for consanguineous
    public boolean RS_ONLY = false; // only consider snps with rsIDs to be snps
    public int REF_CUTOFF = 1;
    public int NUM_EXON_READ = -1;
    public String BAM_FILE_NAMES;
    // TODO: remove the following
    //public String BAM_EXT = ".bam";

    // Mutation Types
    public boolean FRAMESHIFT = false;
    public boolean STOP_GAINED = false;
    public boolean SPLICE_SITE = false;
    public boolean MISSENSE = false;

    // Filters
    public boolean FILTER_QUALITY = true;
    public boolean FILTER_INHERITANCE = true;
    public boolean FILTER_SNP = true;
    public boolean FILTER_MUT_TYPE = true;

    // Output Specifications
    public boolean PRINT_BY_FAM = false;
    public boolean PRINT_SHORT_GENE = false;


    public Flags() {
    }


    /* Sets the Header flags according to header column given.
     * Resets all header flags to -1 before setting again.
     * Meant to be called before each file processed in case files have different
     * headers / columns available.
     */
    public void parseHeader(String[] headerCols) {
        // Reset Header Num Flags
        POS = -1;
        DEPTH = -1;
        PASS = -1;
        SNP = -1;
        GENE = -1;
        HOMO = -1;
        MUT_TYPE = -1;
        TRANSCRIPT = -1;
        GMAF = -1;
        REF = -1;

        for (int i = 0; i < headerCols.length; i++) {
            if (headerCols[i].contains("epth")) {
                DEPTH = i;
            }
            else if (headerCols[i].contains("osition")) {
                POS = i;
            }
            else if (headerCols[i].contains("egion") || headerCols[i].contains("ffect")) {
                MUT_TYPE = i;
            }
            else if (headerCols[i].contains("gene")) {
                GENE = i;
            }
            else if (headerCols[i].contains("type")) {
                HOMO = i;
            }
            else if (headerCols[i].contains("ilter")) {
                PASS = i;
            }
            else if (headerCols[i].contains("dbsnp") || headerCols[i].contains("external")) {
                SNP = i;
            }
            else if (headerCols[i].contains("ranscript")) {
                TRANSCRIPT = i;
            }
            else if (headerCols[i].contains("GMAF")) {
                GMAF = i;
            }
            else if (headerCols[i].equals("Change") || headerCols[i].equals("change")) {
                REF = i;
            }
        }
    }

    /* Parses Commandline arguments (only called once).
     * And sets appropriate instances of Flags.
     */
    public void parseCommandLineArgs(String[] args) {
        // Assume args array contains all of args where first two are input file with individual names and second is output root
        for (int i = 2; i < args.length; i++) {
            String s = args[i];
            // If no mutation type filters are specified, the default is to look for all non-intronic / non-synonymous
            // -Ext arg filters only for frameshift, stop_gained, and splice_site mutations
            if (s.equals("-Ext")) {
                FRAMESHIFT = true;
                STOP_GAINED = true;
                SPLICE_SITE = true;
            }

            // Filter for specific combination of mutations (comma separated list)
            // m = missense, fs = frameshift, ss = splice site, s = stop
            else if (s.equals("-Mut_Type")) {
                if ((i+1) < args.length) {
                    i++;
                    String[] mutTypes = args[i].split(",");
                    for (String m : mutTypes) {
                        if (m.equals("m")) {
                            MISSENSE = true;
                        }
                        else if (m.equals("fs")) {
                            FRAMESHIFT = true;
                        }
                        else if (m.equals("ss")) {
                            SPLICE_SITE = true;
                        }
                        else if (m.equals("s")) {
                            STOP_GAINED = true;
                        }
                    }
                }
                // Write error message / warning message to stdout
                else {
                    System.out.println("WARNING: no mutation type specified for filtering.  Default will be used.");
                }
            }

            // Pass in Depth Quality
            // -d=int
            else if (s.contains("-d=")) {
                // parse integer
                DEPTH_CUTOFF = Integer.parseInt(s.substring(3));
            }

            // X Linked
            // -x
            // Takes precedence over other inheritance patterns so only returns X mutations
            else if (s.equals("-x")) {
                X_LINKED = true;
                RECESSIVE = false;
            }

            // Include homozygous mutations, even if looking at recessive, non-consanguineous or dominant disease
            else if (s.equals("-homo")) {
                INCLUDE_HOMO = true;
            }

            // Recessive Disease
            // -r
            else if (s.equals("-r")) {
                RECESSIVE = true;
            }

            // SNP frequency to filter off of (default is 0)
            // -snpF=#
            else if (s.contains("-snpF=")) {
                SNP_FREQ = Double.parseDouble(s.substring(6));
            }

            // Reference Directory to do Cross Reference off of
            // Assume list of semi-colon separated directories (whole paths)
            // -ref "PATH!;PATH2;..."
            // NOTE: might need extra slash on end of the last filename
            else if (s.equals("-ref")) {
                if ((i+1) < args.length) {
                    i++;
                    String[] dirs = args[i].split(";");
                    CROSS_REF_DIR = new String[dirs.length];
                    for (int j = 0; j < dirs.length; j++) {
                        CROSS_REF_DIR[j] = dirs[j];
                    }
                }
                else {
                    System.out.println("WARNING: no reference directories specified.  No Cross Reference will be Executed.");
                }
            }

            else if (s.equals("-ref_cutoff=")) {
                REF_CUTOFF = Integer.parseInt(s);
            }

            // Consanguinity filter
            // -c=#
            // if # is 0, then filter for all HOMO mutations, but otherwise, filter for HOMO mutations
            // within HOMO Regions of at least # Cm
            else if (s.contains("-c=")) {
                CONSANGUINEOUS = Integer.parseInt(s.substring(3));
            }

            else if (s.equals("-cFam")) {
                if ((i+1) < args.length) {
                    i++;
                    String[] dirs = args[i].split(",");
                    CONSANGUINEOUS_FAMILIES = new ArrayList<Integer>();
                    for (String fam : dirs) {
                        CONSANGUINEOUS_FAMILIES.add(Integer.parseInt(fam));
                    }
                }
                else {
                    System.out.println("WARNING: consanguinity of families not specified.  Consanguineous flag will be applied to all families.");
                }
            }

            // rsIDs used to mark snps only
            // -rs
            else if (s.equals("-rs")) {
                RS_ONLY = true;
            }

            // exon Reads incorporation
            // if none given, assume all individuals must possess the mutation (i.e. no exonRead incorporation)
            // Note: is implemented on a per family basis (i.e. number of individuals per family that must have mutation; same for all families)
            else if (s.contains("-exonRead=")) {
                NUM_EXON_READ = Integer.parseInt(s.substring(10));
            }

            // bam File Directory; print warning if no bam file directory given, but exonRead specified
            else if (s.equals("-bamRef")) {
                if ((i+1) < args.length) {
                    i++;
                    BAM_FILE_NAMES = args[i];
                }
            }

            // bam extension
            /*
            else if (s.equals("-bamExt")) {
                if ((i+1) < args.length) {
                    i++;
                    BAM_EXT = args[i];
                }
            }
            */

            // List of Filters to Use
            else if (s.equals("-filters")) {
                FILTER_INHERITANCE = false;
                FILTER_MUT_TYPE = false;
                FILTER_QUALITY = false;
                FILTER_SNP = false;

                if ((i+1) < args.length) {
                    i++;
                    String[] filters = args[i].split(",");
                    for (String m : filters) {
                        if (m.equals("q")) {
                            FILTER_QUALITY = true;
                        }
                        else if (m.equals("snp")) {
                            FILTER_SNP = true;
                        }
                        else if (m.equals("in")) {
                            FILTER_INHERITANCE = true;
                        }
                        else if (m.equals("mt")) {
                            FILTER_MUT_TYPE = true;
                        }
                    }

                }
            }

            // Output Flag
            else if (s.equals("-outputFam"))
                PRINT_BY_FAM = true;

            else if (s.equals("-outputShortGene"))
                PRINT_SHORT_GENE = true;

            else if (s.equals("-help")) {
                System.out.println("Usage:");
                System.out.println("java FilterProgram fileWithExomeFileNames outputFileRootName [opts]");
                System.out.println();
                System.out.println("Options");
                System.out.println("----------");
                System.out.println("Flag Format \t Description");
                System.out.println("-d= \t Retain mutations sequenced with coverage greater than or equal to d.  Default is 20 for consanguineous and 10 otherwise.");
                System.out.println("-Ext \t Retain mutations of type Frameshift, Splice_Site, or Stop_gained/lost only.");
                System.out.println("-Mut_Type {m,s,fs,ss} \t Retain mutations only of types specified (missense, stop, frameshift, splice_site)");
                System.out.println("-x \t X linked disease. Default is false.  Returns homozygous mutations on X chromosome only.  All input files assumed to be male.");
                System.out.println("-r \t Recessive disease.  Default is Dominant.");
                System.out.println("-c= \t Consanguineous Individuals.  Retain mutations within homozygous regions of specified cM. (Note: applied to all families unless cFam is specified.)");
                System.out.println("-cFam {in1,in2,in3....(in : 0 to #families-1) \t Indicates which families are consanguineous if -c= flag is specified.  Note: same cM distance applied to all consanguineous families.");
                System.out.println("-rs \t Consider only snps with an rsID in filtering for known snps.  Default is to consider all Ids.");
                System.out.println("-snpF= \t Retain mutations, if common snps, with a frequency less than this.  Default is 0.");
                System.out.println("-ref {dir1;dir2;dir3...} \t Reference Directories to cross reference mutations in candidate individuals against.");
                System.out.println("-ref_cutoff= \t Cutoff for Cross Reference.  Eliminate mutations if found in at least this many individuals.  Default is 1.");
                System.out.println("-filters {q,snp,in,mt} \t Specify which filters to use.  Default is to use all filters (quality, snps, inheritance pattern (homo/heterozygous), mutation type (all non-intronic and non-synonymous))");
                System.out.println("-exonRead=# \t Consider mutations that occur in at least # individuals and check bam files to see if these locations are unread in other individuals.");
                System.out.println("-bamRef {filename} \t File name with list of bam files in same format as list of tab-separated files.");
                //System.out.println("-bamRef {dir} \t Read Bam Files from within this directory.  Names of Bam files assumed to be same as roots of *.annot.tab");
                //System.out.println("-bamExt {ext} \t Extension of bam files (following NAME in NAME.annot.tab).  Default is .bam.  Specify root begining with '.'.");
                System.out.println("-outputFam \t Prints potential mutations grouped by family first over grouping by gene.  Default is to group overall by gene.");
                System.out.println("-outputShortGene \t Prints only those genes that more than one family have a mutation in.  If only one family provided, nothing is printed.  Flag ignored if output grouped by family.");
            }
            else {
                System.out.println("WARNING: Flag not Recognized.  Please type -help for further information.");
            }
        }

        // check for Exon Read but no Bam File Directory
        if (NUM_EXON_READ > 0 && (BAM_FILE_NAMES == null)) {
            System.out.println("WARNING: No bam file names specified.  Exon Read cannot be executed.");
        }
    }



}
