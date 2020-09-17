package dataprep;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Provides methods for converting cleaned FAERS data into an arff file
 * @author zhengc
 *
 */
public class WidthToArff {

	public static void main(String[] args) throws IOException {
//		String widthfile = args[0];
//		String arfffile = args[1];

		String widthfile = "/Users/zhengc/workspace/FARES/data/FARES/UMLS_map_data/ID_indications_all_clean_width_umls_id_diso_AD_filtered";
		String arfffile = "/Users/zhengc/workspace/FARES/data/FARES/UMLS_map_data/ID_indications_all_clean_width_umls_id_diso_AD_filtered_sp.arff";
				
		// write header info
		BufferedReader br = new BufferedReader(new FileReader(new File(widthfile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(arfffile)));
		bw.write("@relation " + "widthfile" + "\n");
		bw.write("\n");
		
//		Map<String, Integer> itemindex = getItemIndex(csvfile);
		Set<String> itemset = getItemSet(widthfile);
		System.out.println(itemset.size());
		List<String> itemlist = new ArrayList<String>(itemset);
		Collections.sort(itemlist);
		for (String item : itemlist) {
			bw.write("@attribute " + "'" + item + "'" + " {0,1}" + "\n");
		}
		
		// write data
		bw.write("\n");
		bw.write("@data" + "\n");
//		String line = null;
		String line = br.readLine();
//		while ((line=br.readLine()) != null) {
//			bw.write(line + "\n");
//		}
		
		while ((line=br.readLine()) != null) {
			StringBuilder sb = new StringBuilder();
			String[] parts = line.split(",");
			
			// get unique indications
			Set<String> partsset = new HashSet<String>();
			for (int i=0; i<parts.length; i++) {
				partsset.add(parts[i]);
			}
			// convert to list and sort the list
			List<String> partslist = new ArrayList<String>(partsset);
			Collections.sort(partslist);
			
			bw.write("{");
			for (int i=0; i<partslist.size(); i++) {
				int index = itemlist.indexOf(partslist.get(i));
				sb.append(index + " " + "1" + ",");
			}
			sb.append("}");
			String spline = sb.toString(); 
			String splinef = spline.replace(",}", "}"); //remove last comma
			bw.write(splinef + "\n");
		}
		bw.close();
		br.close();
	}
	
	private static Map<String, Integer> getItemIndex(String filename) throws IOException {
		Map<String,Integer> itemindex = new HashMap<String, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		int m = 0;
		while ((line=br.readLine()) != null) {
			String[] parts = line.split(",");
			for (int i = 0; i < parts.length; i++) {
				if (!itemindex.containsKey(parts[i])) {
					itemindex.put(parts[i], m);
					m++;
				}
			}
		}
		br.close();
		return itemindex;
	}
	
	private static Set<String> getItemSet(String filename) throws IOException {
		Set<String> itemset = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
//		String line = null;
		String line = br.readLine();
		while ((line=br.readLine()) != null) {
			String[] parts = line.split(",");
			for (int i = 0; i < parts.length; i++) {
					itemset.add(parts[i]);
					
			}
		}
		br.close();
		return itemset;
	}
}

