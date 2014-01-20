package com.example.mutationfilter;

import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/10/13
 * Time: 2:22 PM
 * To change this template use File | Settings | File Templates.
 *
 * Filters Mutations based on Homo/Heterozygosity, based on the assumed inheritance pattern as given.
 * If Family is Consanguineous, will look for sections of homozygosity of at least a certain length in cM.
 *
 * Note: Homo/Heterozygosity is ignored with X chromosome since calls are generally incorrect with this.
 *       Technically, mutation in X, if x-linked disease, will be HET since male only has one copy of X.
 *
 */
public class InheritanceFilter implements Filter {

    // Constants
    private final double PERCENT = 1; // percentage cutoff for homozygous regions

    // Instance Variables
    private Flags flags;
    private HomozygousRegions homoRegions; // set this separately for each person; only used for Consanguineous Filtering (> 0)
    private CmToMbMap conversion; // set once and used for all individuals; only used for Consanguineous Filtering (> 0)

    public InheritanceFilter(Flags f) throws IOException {
        flags = f;
        homoRegions = null;
        // assume if consanguineous flag is set that at least one family is consanguineous
        if (f.CONSANGUINEOUS > 0) {
            conversion = new CmToMbMap();
        }
        else {
            conversion = null;
        }
    }

    // Allows the flag variable to be reset each time another individual is seen
    // Allows for flexibility in information that each person might have available to filter on
    public void updateFlagHeaders(String[] header) {
        flags.parseHeader(header);
    }

    public boolean pass(String[] tabRow, FamilyDataGroup family) {
        // if HOMO header not set or not found, location cannot be filtered out
        if (flags.HOMO < 0) {
            return true;
        }

        // only keep X locs if X_LINKED is true; otherwise, keep no X locs
        // TODO: POSSIBLY CHANGE THIS IF GET TOO MANY ADDITIONAL LOCATIONS
        boolean isX = tabRow[0].split(":")[0].equalsIgnoreCase("x");
        if (flags.X_LINKED) {
            if (isX && tabRow[flags.HOMO].contains("HOMO")) {
                return true;
            }
            return false;
        }

        else {
            if (isX) {
                return false;
            }
        // if dominant or recessive, but non-consanguineous, look for HET mutations only
        if (!flags.RECESSIVE || !(family.consanguineous)) {
            if (!tabRow[flags.HOMO].contains("HOMO")) {
                return true;
            }
            return false;
        }
        else {
            if (flags.CONSANGUINEOUS == 0) {
                if (tabRow[flags.HOMO].contains("HOMO")) {
                    return true;
                }
                return false;
            }
            // Filter by homozygous regions rather than individually homozygous locations
            else {
                if (conversion == null || homoRegions == null) {
                    System.err.println("ERROR: cMToMb map or homozygous region for specific individual not set.  Consanguineous Filtering Could not be Completed.");
                    return false;
                }
                // Assume the homozygous Region object is set to correct person
                String[] locParts = tabRow[0].split(":");
                String c = locParts[0];
                int loc = Integer.parseInt(locParts[1]);
                if (homoRegions.containsLoc(c, loc)) {
                    return true;
                }
                return false;
            }

        }
        }
    }

    // Set Homozygous Regions for specific person
    public void setHomoRegions(String filename) throws IOException {
        homoRegions = new HomozygousRegions(PERCENT, flags.CONSANGUINEOUS, flags, filename, conversion);
    }
}
