package com.example.mutationfilter;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/6/13
 * Time: 3:17 PM
 * To change this template use File | Settings | File Templates.
 *
 * Filters Mutations based on their type (consequence/effect).  Default is to remove all mutations
 * that are intronic, intergenic, synonymous
 */
public class MutationTypeFilter implements Filter {

    // Instance Variables
    private Flags flags;

    public MutationTypeFilter(Flags flags) {
        this.flags = flags;
    }

    // Allows the flag variable to be reset each time another individual is seen
    // Allows for flexibility in information that each person might have available to filter on
    public void updateFlagHeaders(String[] header) {
        flags.parseHeader(header);
    }

    public boolean pass(String[] tabRow, FamilyDataGroup family) {
        if (flags.MUT_TYPE < 0) {
            return true;
        }

        String intron = "intron";
        String downstream = "downstream";
        String upstream = "upstream";
        String nonCoding = "non_coding";
        String interGene = "intergen";
        String synonymous = "synonymous";
        String nmd_transcript = "nmd_transcript";
        String nc_transcript = "nc_transcript";

        String actualMut = tabRow[flags.MUT_TYPE].toLowerCase();

        // Default (none of mut_type commandline args are set)
        if (!flags.FRAMESHIFT && !flags.STOP_GAINED && !flags.SPLICE_SITE && !flags.MISSENSE) {
            if ((actualMut.contains(intron) && !actualMut.contains("splice")) || actualMut.contains(downstream)
                    || actualMut.contains(upstream) || actualMut.contains(nonCoding)
                    || actualMut.contains(interGene) || actualMut.contains(synonymous)
                    || actualMut.contains(nmd_transcript) || actualMut.contains(nc_transcript)) {
                return false;
            }
            return true;
        }

        // Specified type of mutation
        boolean useFrameShift = (flags.FRAMESHIFT && actualMut.contains("frameshift"));
        boolean useStop = (flags.STOP_GAINED && actualMut.contains("stop"));
        boolean useSplice = (flags.SPLICE_SITE && actualMut.contains("splice"));
        boolean useMissense = (flags.MISSENSE && actualMut.contains("missense"));

        if (useFrameShift || useStop || useSplice || useMissense) {
            return true;
        }
        else {
            return false;
        }
    }
}
