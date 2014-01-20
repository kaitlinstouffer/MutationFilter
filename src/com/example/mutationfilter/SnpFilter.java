package com.example.mutationfilter;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 3:50 PM
 * To change this template use File | Settings | File Templates.
 *
 * Filters mutations based on their being registered snps and/or their frequency if
 * they're registered snps.
 */
public class SnpFilter implements Filter {

    // Instance Variables
    private Flags flags;   // Used for header Nums

    public SnpFilter(Flags flags) {
        this.flags = flags;
    }

    // Allows the flag variable to be reset each time another individual is seen
    // Allows for flexibility in information that each person might have available to filter on
    public void updateFlagHeaders(String[] header) {
        flags.parseHeader(header);
    }

    public boolean pass(String[] tabRow, FamilyDataGroup family) {
        // don't filter based off snp if column not there
        if (flags.SNP < 0) {
            return true;
        }
        double f = -1;
        boolean isSNP = false;
        boolean freqSet = true;

        // Considers rs and any other IDs as SNPs
        if (flags.RS_ONLY) {
            if (((flags.SNP < tabRow.length && tabRow[flags.SNP] != null) && !(tabRow[flags.SNP]).isEmpty()) &&
                    ((tabRow[flags.SNP].contains("rs")))) {
                f = 0;
                isSNP = true;
                freqSet = false;
            }

        }
        else {
            if (((flags.SNP < tabRow.length && tabRow[flags.SNP] != null) && !(tabRow[flags.SNP]).isEmpty()) &&
                    ((tabRow[flags.SNP].contains("rs")) || (!tabRow[flags.SNP].equals(".")))) {
                f = 0;
                isSNP = true;
                freqSet = false;
            }
        }

        if (((flags.GMAF > 0) && flags.GMAF < tabRow.length) && !(tabRow[flags.GMAF] == null || tabRow[flags.GMAF].isEmpty())) {
            String allele = null;
            if (flags.REF > 0) {
                allele = tabRow[flags.REF].split(">")[1];
            }
            String[] alls = tabRow[flags.GMAF].split("&");
            for (String s : alls) {
                if (s.startsWith(allele)) {
                    f = Double.parseDouble(s.split(":")[1]);
                    freqSet = true;
                    break;
                }
            }

            // if frequency is not set, then subtract all other frequencies from 1
            if (isSNP && !freqSet) {
                f = 1;
                for (String s : alls) {
                    f -= Double.parseDouble(s.split(":")[1]);
                }
            }

        }

        // Compare to command line frequency
        if (f < flags.SNP_FREQ) {
            return true;
        }
        return false;
    }
}
