package com.example.mutationfilter;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: kaitlinstouffer
 * Date: 12/5/13
 * Time: 11:18 AM
 * To change this template use File | Settings | File Templates.
 *
 * Class with methods to facilitate reading from tab-separated file.
 */
public class TabFileReader extends BufferedReader {
    // Inherits all instance variables and constructors
    public TabFileReader(String fileName) throws FileNotFoundException {
        super(new FileReader(fileName));
    }

    // Reads in Line and returns list of Elements
    public List<String> readColumns() throws IOException {
        String line = super.readLine();
        if (line == null)
            return null;
        else {
            List<String> columns = Arrays.asList(line.split("\t"));
            return columns;
        }
    }

    // Reads in Line and returns array of Elements
    public String[] readColumnsArray() throws IOException {
        String line = super.readLine();
        if (line == null)
            return null;
        else {
            return line.split("\t");
        }
    }
}
