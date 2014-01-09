package com.example.mutationfilter;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/6/13
 * Time: 3:04 PM
 * To change this template use File | Settings | File Templates.
 *
 * Filters Mutations by their DEPTH and a PASSING filter (only elminates HARD_TO_VALIDATE results)
 */
public class QualityFilter implements Filter {

    // Instance Variables
    private Flags flags;

    public QualityFilter(Flags flags) {
        this.flags = flags;
    }

    // Allows the flag variable to be reset each time another individual is seen
    // Allows for flexibility in information that each person might have available to filter on
    public void updateFlagHeaders(String[] header) {
        flags.parseHeader(header);
    }

    public boolean pass(String[] tabRow) {
        if (flags.DEPTH < 0) {
            return true;
        }

        int depthCutOff = flags.DEPTH_CUTOFF;

        // Filter off of both Depth and PASSING filter (depth must be greater than cutoff and PASS must not be hard to validate)
        if (flags.DEPTH_CUTOFF < 0) {
            if (flags.CONSANGUINEOUS >= 0) {
                depthCutOff = 20;
            }
            else {
                depthCutOff = 10;
            }
        }

        if (Integer.parseInt(tabRow[flags.DEPTH]) >= depthCutOff) {
            if (flags.PASS > 0) {
                if (!tabRow[flags.PASS].equalsIgnoreCase("HARD_TO_VALIDATE")) {
                     return true;
                }
                return false;
            }
            return true;
        }
        return false;

    }
}
