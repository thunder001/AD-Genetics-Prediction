package dataprep;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**  
 * Provides methods for processing UMLS mapped FAERS data
 * @author zhengc
 *
 */
public class ProcessUMLS {

	public static void main(String[] args) throws IOException {
		String input = "/Users/zhengc/workspace/FARES/data/FARES/UMLS_map_data/ID_indications_all_clean_width_umls_id_diso_AD";
		String output = "/Users/zhengc/workspace/FARES/data/FARES/UMLS_map_data/ID_indications_all_clean_width_umls_id_diso_AD_filtered";
		String idnummapfile = "/Users/zhengc/workspace/FARES/data/FARES/UMLS_map_data/ID_indications_all_clean_width_umls_id_diso_AD_idnummap";
		
		// Compute total terms
		Set<String> terms = getTermSet(input);
		System.out.println(terms.size());
		// Create numbers of each termID appeared in the data
		Map<String, Integer> idnummap = getIdNumMap(input);
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(idnummapfile)));
		bw.write("UMLS_id" + "," + "Number" + "\n");
		for (String id : idnummap.keySet()) {
			bw.write(id + "," + idnummap.get(id) + "\n");
		}
		bw.close();
		
		System.out.println("Number of total diseases: " + idnummap.size());
		
		// Set threshold and get processed data
		Set<String> freqset = getFreqTermSet(idnummap, 1, 500000);
		System.out.printf("%d disease will be filtered", freqset.size());
//		for (String term : freqset) {
//			System.out.println(term);
//		}
		procUMLS(input, freqset, output);
	}
	
	/**
	 * Get all disease terms
	 * @param filename raw FDA labeling data file
	 * @return disease term set
	 * @throws IOException
	 */
	private static Set<String> getTermSet(String filename) throws IOException {
		Set<String> termSet = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line=br.readLine()) != null) {
			String[] terms = line.split(",");
			for (int i=0; i<terms.length; i++) {
				termSet.add(terms[i]);
			}
		}
		br.close();
		return termSet;
	}
	
	/**
	 * Compute the number of occurrence of each disease term in FAERS data 
	 * @param filename cleaned, UMLS mapped FAERS file
	 * @return id_number map that key is disease UMLS term and value is number of occurrence
	 * @throws IOException
	 */
	private static Map<String, Integer> getIdNumMap(String filename) throws IOException {
		Map<String, Integer> idNumMap = new HashMap<String, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line=br.readLine()) != null) {
			String[] terms = line.split(",");
			// get unique terms
			Set<String> termset = new HashSet<String>();
			for (int i=0; i<terms.length; i++) {
				termset.add(terms[i]);
			}
			for (String term : termset) {
				if (!idNumMap.containsKey(term)) {
					idNumMap.put(term, 1);
				} else {
					idNumMap.put(term, idNumMap.get(term) + 1);
					
				}
			}
		}
		br.close();
		return idNumMap;
	}
	
	/**
	 * Remove rare and too many occurrence disease term in FAERS
	 * @param idNumMap id_number map
	 * @param minfreq low bound of threshold
	 * @param maxfreq up bound threshold
	 * @return filtered disease term set 
	 */
	private static Set<String> getFreqTermSet(Map<String, Integer> idNumMap, int minfreq, int maxfreq) {
		Set<String> freqSet = new HashSet<String>();
		int value = 0;
		for (String key : idNumMap.keySet()) {
			value = idNumMap.get(key);
			if (value <= minfreq | value >= maxfreq) {
				freqSet.add(key);
			}
		}
		return freqSet;
	}
	
	/**
	 * Gets filtered FAERS data
	 * @param input cleaned, UMLS mapped FARES data
	 * @param freqset filtered disease term set 
	 * @param output filter FAERS data
	 * @throws IOException
	 */
	private static void procUMLS(String input, Set<String> freqset,  String output) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(input)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output)));
		String line = null;
		while ((line=br.readLine()) != null) {
			String[] terms = line.split(",");
			// get unique terms
			Set<String> termset = new HashSet<String>();
			for (int i=0; i<terms.length; i++) {
				termset.add(terms[i]);
			}
			List<String> termlist = new ArrayList<String>(termset);
			
			StringBuilder sb = new StringBuilder();
			for (String term : termlist) {
				if (!freqset.contains(term)) {
					sb.append(term + ",");
				}
			}
			String sbline = sb.toString();
//			System.out.println(sbline);
			if (!sbline.isEmpty()) {
				sbline = sbline.replaceFirst(",$", "");
				bw.write(sbline + "\n");
			}
		}
		bw.close();
		br.close();
	}

}
