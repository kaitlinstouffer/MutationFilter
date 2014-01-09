package com.example.mutationfilter;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 11:55 AM
 * To change this template use File | Settings | File Templates.
 *
 * Class for facilitating writing to tab separated files.
 */
public class TabFileWriter extends BufferedWriter {

    // Inherits constructor and Instance Variables
    public TabFileWriter(String fileName) throws IOException {
        super(new FileWriter(new File(fileName)));
    }

    // Writes list of Elements
    public void writeColumns(List<String> cols) throws IOException {
        for (String s : cols) {
            super.write(s + "\t");
            super.write("\n");

        }
    }

    // Writes array of elements
    public void writeColumns(String[] cols) throws IOException {
        for (String s : cols) {
            super.write(s + "\t");
            super.write("\n");
        }
    }


}
