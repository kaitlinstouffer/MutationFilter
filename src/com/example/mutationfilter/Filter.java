package com.example.mutationfilter;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 3:42 PM
 * To change this template use File | Settings | File Templates.
 *
 * Filter Interface to define functions of general filter.
 * Implemented by: SnpFilter.java, HomozygosityFilter.java, MutationTypeFilter.java, QualityFilter.java
 */
public interface Filter {

    /* Returns true if Mutation specified by input row of .tab file passes the specifications
     * of the given filter.  FamilyDataGroup family indicates the FamilyDataGroup to which the individual
     * whose row is being filtered belongs to.
     *
     * Returns false if position does not pass.
     */
    boolean pass(String[] tabRow, FamilyDataGroup family);

    /* Updates the Flags variable assumed to be in each implementation of Filter by
     * setting the header to whatever is given by
     * String[] header.
     */
    void updateFlagHeaders(String[] header);


}
