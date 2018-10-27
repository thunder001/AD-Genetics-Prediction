package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Utility class to read a txt file into a List
 * @author zhengc
 *
 */
public class FileToList {

	/**
	 * Creates a List from a txt file
	 * @param filename file that contains items
	 * @return a List
	 * @throws IOException
	 */
	public static List<String> createListFromfile(String filename) throws IOException {
		List<String> genelist = new ArrayList<>();
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String gene = null;
		while ((gene = br.readLine()) != null) {
			genelist.add(gene);
		}
		br.close();
		return genelist;
	}

}
