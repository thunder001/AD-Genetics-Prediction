package util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * This class provides disease mapping between DCN and OMIM disease
 * @author zhengc
 *
 */
public class DCNOMIMUMLSIDmap {

	public static Map<String, List<String>> dcnidnamemap = new HashMap<>();
	public static Map<String, String> dcnnameidmap = new HashMap<>();
	public static Map<String, List<String>> omimidnamemap = new HashMap<>();
	public static Map<String, String> omimnameidmap = new HashMap<>();
	
	private static Map<String, List<String>> DCN_OMIM = new HashMap<>();
	
	/**
	 * Initiates and creates a DCNOMIMUMLSIDmap object
	 * @param DCN_disease DCN disease name
	 * @param OMIM_dispattern OMIM disease name pattern that equals DCN disease name
	 */
	public DCNOMIMUMLSIDmap(String DCN_disease, String OMIM_dispattern) {
		
		String DCN_id = dcnnameidmap.get(DCN_disease);
		
		List<String> OMIM_ids = new ArrayList<>();
		for (String dis : omimnameidmap.keySet()) {
			if (dis.contains(OMIM_dispattern)) {
				OMIM_ids.add(omimnameidmap.get(dis));
			}	
		}
		Collections.sort(OMIM_ids);
		DCN_OMIM.put(DCN_id, OMIM_ids);
		
	}
	
	/**
	 * Creates a DCN name_UMLS mapping
	 * @param mapfile a name_UMLS mapping file
	 * @throws IOException
	 */
	public static void createDCNIdNameMap(String mapfile) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(mapfile));
		String line = null;
		while((line = br.readLine()) != null){
			String[] parts = line.split("\\|");
			String name = parts[0].toLowerCase();
			String id = parts[1];
			if (!dcnnameidmap.containsKey(name)) {
				dcnnameidmap.put(name, id);
			}
			
			if (!dcnidnamemap.containsKey(id)) {
				List<String> names = new ArrayList<>();
				names.add(name);
				dcnidnamemap.put(id, names);
			} else {
				dcnidnamemap.get(id).add(name);
			}
		}
		br.close();
	} 
	

	/**
	 * Creates a OMIM name_UMLS mapping
	 * @param mapfile a name_UMLS mapping file
	 * @throws IOException
	 */
	public static void createOMIMIdNameMap(String mapfile) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(mapfile));
		String line = null;
		while((line = br.readLine()) != null){
			String[] parts = line.split("\\|");
			String name = parts[0].toLowerCase();
			String id = parts[1];
			if (!omimnameidmap.containsKey(name)) {
				omimnameidmap.put(name, id);
			}
			
			if (!omimidnamemap.containsKey(id)) {
				List<String> names = new ArrayList<>();
				names.add(name);
				omimidnamemap.put(id, names);
			} else {
				omimidnamemap.get(id).add(name);
			}
		}
		br.close();
	} 
	
	/**
	 * Get DCN_OMIM mapping
	 * @return DCN_OMIM mapping
	 */
	public Map<String, List<String>> getDCNOMIMUMLSIDmap() {
		return DCN_OMIM;
	}
	
	public static void main(String[] args) throws IOException {

		String dcnmapfile = "/Users/zhengc/workspace/FARES_final/data/FARES/map/term_umls_id_diso";
		String omimmapfile = "/Users/zhengc/workspace/FARES_final/data/OMIM/map/OMIM_umls_id_diso";
		createDCNIdNameMap(dcnmapfile);
		createOMIMIdNameMap(omimmapfile);
		
		String DCN_dis = "dementia";
		String OMIM_pat = "alzheimer";
		
		DCNOMIMUMLSIDmap domim = new DCNOMIMUMLSIDmap(DCN_dis, OMIM_pat);
		Map<String, List<String>> dcn_omim = domim.getDCNOMIMUMLSIDmap();
				
		System.out.println(DCN_dis + ": " + dcnnameidmap.get(DCN_dis));
		
		String DCN_dis_id = dcnnameidmap.get(DCN_dis);
		List<String> OMIM_ids = dcn_omim.get(DCN_dis_id);
		
		for (String id : OMIM_ids) {
			List<String> names = omimidnamemap.get(id);
			for (String name : names) {
				System.out.println(name + ": " + id);
			}
		}

	}

}
