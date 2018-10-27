package network;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import GraphAlgorithm.ComSearch;
import util.CollectionsEx;
import util.DCNOMIMUMLSIDmap;

/**
 * Provides methods for creating a bipartite network and write it to a network file
 * @author zhengc
 *
 */
public class CommGeneticsLoader {
	public static HashSet<String> entry_set = new HashSet<String>();
	public static HashSet<String> disease_set = new HashSet<String>();
	public static HashSet<String> gene_set = new HashSet<String>();
	public static HashSet<String> DMN_set = new HashSet<String>();
	public static HashSet<String> DCN_set = new HashSet<String>();
	
	public static List<String> entry_list = new ArrayList<String>();
	public static List<String> disease_list = new ArrayList<String>();
	public static List<String> gene_list = new ArrayList<String>();
	public static List<String> DMN_list = new ArrayList<String>();
	public static List<String> DCN_list = new ArrayList<String>();
	
	public static int num_gene = 0;
	public static int num_dmn;
	public static int num_dcn;
	public static int num_disease;
	
	public static Map<Integer, Integer> dis_gene = new HashMap<Integer, Integer>();
	
	
	public static Map<String, Integer> entry_index = new HashMap<String, Integer>();
	public static Map<Integer, String> index_entry = new HashMap<Integer, String>();
	
	public static Map<String, String> idnamemap = new HashMap<String, String>();
	public static Map<String, String> nameidmap = new HashMap<String, String>();
	public static Map<String, String> OMIMidnamemap = new HashMap<String, String>();
	public static Map<String, String> OMIMnameidmap = new HashMap<String, String>();
	
	public static HashSet<String> soc_set = new HashSet<String>();
	public static Map<String, String> umlssocmap = new HashMap<String, String>();
	public static Map<String, Integer> soc_idx = new HashMap<String, Integer>();
	public static Map<Integer, String> idx_soc = new HashMap<Integer,String>();
	public static final double P = 0.85;
	
	/**
	 * Read DCN net and get DCN disease list
	 * @param filename DCN net file
	 * @throws IOException
	 */
	public static void readCommNet(String filename) throws IOException{
		// get all unique disease set
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			String disease = line.split("\\|")[0];
			DCN_set.add(disease);
		}
		br.close();
		DCN_list = util.CollectionsEx.setToList(DCN_set);
		num_dcn = DCN_set.size();
		System.out.println("number of DCN: " + num_dcn);
		
	}
	
	/**
	 * Write a map from DCN disease UMLS to disease name
	 * @param imput map file from disease UMLS to name 
	 * @param output a file to be written
	 * @throws IOException
	 */
	public static void writeDCNIdName(String imput, String output) throws IOException {
		Map<String, String> idnamemap = new HashMap<String, String>();
		
		BufferedReader br = new BufferedReader(new FileReader(new File(imput)));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			String disease_id = line.split("\\|")[0];
			String disease_name = line.split("\\|")[1];
			if(!idnamemap.keySet().contains(disease_id)) {
				idnamemap.put(disease_id, disease_name);
			} 
		}
		br.close();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output)));
		Set<String> ids = idnamemap.keySet();
		List<String> id_list = util.CollectionsEx.setToList(ids);
		Collections.sort(id_list);
		
		for(String id : id_list) {
			String name = idnamemap.get(id).replace("-- ", "");
			bw.write(id + "|" + name + "\n");
		}
		bw.close();
	}
	
	/**
	 * Reads DMN to get DMN disease list
	 * @param filename DMN net file
	 * @throws IOException
	 */
	public static void readDMN(String filename) throws IOException{
		// get all unique disease set
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			String disease1 = line.split("\\|")[0];
			DMN_set.add(disease1);
			String disease2 = line.split("\\|")[1];
			DMN_set.add(disease2);
		}
		br.close();
		DMN_list = util.CollectionsEx.setToList(DMN_set);
		num_dmn = DMN_set.size();
		System.out.println("number of DMN: " + num_dmn);
		
	}
	
	/**
	 * Reads PPI to get gene list
	 * @param filename
	 * @throws IOException
	 */
	public static void readPPI(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line=br.readLine()) != null) {
			String gene = line.split("\\|")[0];
//			System.out.println(2);
			gene_set.add(gene);
		}
		br.close();
		gene_list = util.CollectionsEx.setToList(gene_set);
		num_gene = gene_set.size();
		System.out.println("number of gene: " + num_gene);
	}
	
	/**
	 * Reads UMLS-SOC map file to get UMLS-SOC map
	 * @param filename
	 * @throws IOException
	 */
	public static void readUMLSSOC(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line=br.readLine()) != null) {
			String[] parts = line.split("\\$");
			String umls = parts[0];
			String soc = parts[8];
			soc_set.add(soc);
			umlssocmap.put(umls, soc);
		}
		br.close();
		
		int m = 0;
		for (String soc : soc_set) {
			soc_idx.put(soc, m);
			idx_soc.put(m, soc);
			m++;
		}
		
	}
	

	
	
	/**
	 * Creates entry index for all nodes
	 * @param entry_list
	 */
	public static void createIndex(List<String> entry_list) {
		int m = 0;
		for (String e : entry_list) {
			if (!entry_index.containsKey(e)) {
				entry_index.put(e, m);
				index_entry.put(m, e);
				m++;
			}
		}
		System.out.println("Total nodes: " + entry_index.size());
		System.out.println("--------");
	}
	
	/**
	 * Build a bipartite graph including DCN and PPI, additional DCN-OMIM disease mapping also included
	 * @param commnetfile DCN file
	 * @param ppifile PPI file
	 * @param disgenefile disease genetics file from OMIM
	 * @param dcn_omim mapping from DCN disease to OMIM disease
	 * @return a bipartite graph
	 * @throws IOException
	 */
	
	
	public static DisGraph createCommGeneticsGraph(String commnetfile, 
			String ppifile, String disgenefile, Map<String, List<String>> dcn_omim) throws IOException {
		
		readCommNet(commnetfile);
		readPPI(ppifile);
		disease_set.addAll(DCN_set);
		disease_list = util.CollectionsEx.setToList(disease_set);
		entry_list.addAll(disease_list);
		
		entry_list.addAll(gene_list);
		createIndex(entry_list);
		num_disease = disease_list.size();
		
		int n = entry_index.size();
		DisGraph cgg = new DisGraph(n);
		
		
		// construct DCN net from network file
		// Note: there are two type of network files, attention should be paid that one file contains 
		// additional abbreviation of SOC. Consequently, the second entry should be in column 6 instead 5
		BufferedReader br1 = new BufferedReader(new FileReader(new File(commnetfile)));
		String line1 = null;
		while ((line1 = br1.readLine()) != null) {
			String [] s = line1.split("\\|");
//			System.out.println(s[0]);
			//	System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[5])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[5]);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
				
			}
		}
		br1.close();
		
		// Check if diagonal entry is 1
//		for (int j=0; j<disease_list.size(); j++) {
//			System.out.println(cgg.getNet().get(j, j));
//		}
		
		
		BufferedReader br2 = new BufferedReader(new FileReader(new File(ppifile)));
		String line2 = null;
		while ((line2 = br2.readLine()) != null) {
			String [] s = line2.split("\\|");
//			System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[1])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[1]);
//				System.out.println(s[0] +"  " + entry1);
//				System.out.println(s[1] +"  " + entry2);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
			}
		}
		br2.close();
		
		// Check if diagonal entry is 1
//		for (int j=disease_list.size(); j<entry_list.size(); j++) {
//			System.out.println(cgg.getNet().get(j, j));
//		}
		
		BufferedReader br3 = new BufferedReader(new FileReader(new File(disgenefile)));
		String line3 = null;
		Set<String> mappedset = new HashSet<>();
		
		while ((line3 = br3.readLine()) != null) {
			String [] s = line3.split("\\|");
//			System.out.println(s[0]);
			String disid = s[0];
			String gene = s[1];
			if (entry_index.containsKey(disid) && entry_index.containsKey(gene)) {
				mappedset.add(disid);
				int entry1 = entry_index.get(disid);
				int entry2 = entry_index.get(gene);
//					System.out.println(s[0] +"  " + entry1);
//					System.out.println(s[1] +"  " + entry2);
//					double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
			}
			
			// add additional dcn_omim map
			for (String id : dcn_omim.keySet()) {
				List<String> omim_ids = dcn_omim.get(id);
				if (!entry_index.containsKey(disid) && omim_ids.contains(disid)) {
					disid = id; // update id
					if (entry_index.containsKey(disid) && entry_index.containsKey(gene)) {
						int entry1 = entry_index.get(disid);
						int entry2 = entry_index.get(gene);
						cgg.addEdge(entry1, entry2);
					}
				}
			}
			
		}
		br3.close();
		
		// Check if diagonal entry is 1
//		for (int j=disease_list.size(); j<entry_list.size(); j++) {
//			System.out.println(cgg.getNet().get(j, j));
//		}
		
		System.out.println("Total mapped disease: " + mappedset.size());
		
//		for (int i=0; i<n; i++) {
//			cgg.addEdge(i, i);
//		}
		return cgg;
	}
	
	/**
	 * Build a bipartite graph including DCN, DMN and PPI
	 * @param commnetfile DCN file
	 * @param dmnfile DMN file
	 * @param ppifile PPI file
	 * @param disgenefile disease genetics file from OMIM
	 * @return a bipartite graph
	 * @throws IOException
	 */
	public static DisGraph createCommGeneticsGraph(String commnetfile, String dmnfile, 
			String ppifile, String disgenefile) throws IOException {
		
		
		readCommNet(commnetfile);
		readDMN(dmnfile);
		readPPI(ppifile);
		disease_set.addAll(DCN_set);
		disease_set.addAll(DMN_set);
		disease_list = util.CollectionsEx.setToList(disease_set);
		entry_list.addAll(disease_list);
		entry_list.addAll(gene_list);
		createIndex(entry_list);
		num_disease = disease_list.size();
		
		int n = entry_index.size();
		DisGraph cgg = new DisGraph(n); 
//		System.out.println(n);
		
		BufferedReader br1 = new BufferedReader(new FileReader(new File(commnetfile)));
		String line1 = null;
		while ((line1 = br1.readLine()) != null) {
			String [] s = line1.split("\\|");
//			System.out.println(s[0]);
			//	System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[4])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[4]);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry2, entry1);
			}
		}
		br1.close();
		
		BufferedReader br2 = new BufferedReader(new FileReader(new File(dmnfile)));
		String line2 = null;
		while ((line2 = br2.readLine()) != null) {
			String [] s = line2.split("\\|");
//			System.out.println(s[0]);
			//	System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[1])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[1]);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
				
			}
		}
		br2.close();
		
		BufferedReader br3 = new BufferedReader(new FileReader(new File(ppifile)));
		String line3 = null;
		while ((line3 = br3.readLine()) != null) {
			String [] s = line3.split("\\|");
//			System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[1])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[1]);
//				System.out.println(s[0] +"  " + entry1);
//				System.out.println(s[1] +"  " + entry2);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
			}
		}
		br3.close();
		
		BufferedReader br4 = new BufferedReader(new FileReader(new File(disgenefile)));
		String line4 = null;
		Set<String> mappedset = new HashSet<>();
		
		while ((line4 = br4.readLine()) != null) {
			String [] s = line4.split("\\|");
//			System.out.println(s[0]);
			String disid = s[0];
			String gene = s[1];
			if (entry_index.containsKey(disid) && entry_index.containsKey(gene)) {
				mappedset.add(disid);
				int entry1 = entry_index.get(disid);
				int entry2 = entry_index.get(s[1]);
//					System.out.println(s[0] +"  " + entry1);
//					System.out.println(s[1] +"  " + entry2);
//					double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
			}
			
		}
		br4.close();
		
		System.out.println("Total mapped disease: " + mappedset.size());
		
		for (int i=0; i<n; i++) {
			cgg.addEdge(i, i);
		}
		return cgg;
	}
	
	/**
	 * Build a bipartite graph including DCN and PPI, extra mapping and DMN are optional
	 * 
	 * @param commnetfile DCN file
	 * @param dmnfile DMN file
	 * @param ppifile PPI file
	 * @param disgenefile disease genetics file from OMIM
	 * @param dcn_omim mapping from DCN disease to OMIM disease
	 * @param dmn a flag to indicate if DMN will be included
	 * @param extra a flag to indicate if extra mapping will be included
	 * @return a bipartite graph
	 * @throws IOException
	 */
	public static DisGraph createCommGeneticsGraph(String commnetfile, String dmnfile, 
			String ppifile, String disgenefile, Map<String, List<String>> dcn_omim, Boolean dmn, Boolean extra) throws IOException {
		
		
		readCommNet(commnetfile);
		disease_set.addAll(DCN_set);
		if(dmn) {
			readDMN(dmnfile);
			disease_set.addAll(DMN_set);
		}
		disease_list = util.CollectionsEx.setToList(disease_set);
		readPPI(ppifile);
		
		entry_list.addAll(disease_list);
		entry_list.addAll(gene_list);
		createIndex(entry_list);
		num_disease = disease_list.size();
		
		int n = entry_index.size();
		DisGraph cgg = new DisGraph(n); 
//		System.out.println(n);
		
		BufferedReader br1 = new BufferedReader(new FileReader(new File(commnetfile)));
		String line1 = null;
		while ((line1 = br1.readLine()) != null) {
			String [] s = line1.split("\\|");
//			System.out.println(s[0]);
			//	System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[4])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[4]);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry2, entry1);
			}
		}
		br1.close();
		
		if(dmn) {
			BufferedReader br2 = new BufferedReader(new FileReader(new File(dmnfile)));
			String line2 = null;
			while ((line2 = br2.readLine()) != null) {
				String [] s = line2.split("\\|");
//				System.out.println(s[0]);
				//	System.out.println(s[0]);
				if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[1])) {
					int entry1 = entry_index.get(s[0]);
					int entry2 = entry_index.get(s[1]);
//					double score = Integer.parseInt(s[2]);
					cgg.addEdge(entry1, entry2);
					
				}
			}
			br2.close();
		}
		
		
		BufferedReader br3 = new BufferedReader(new FileReader(new File(ppifile)));
		String line3 = null;
		while ((line3 = br3.readLine()) != null) {
			String [] s = line3.split("\\|");
//			System.out.println(s[0]);
			if (entry_index.containsKey(s[0]) && entry_index.containsKey(s[1])) {
				int entry1 = entry_index.get(s[0]);
				int entry2 = entry_index.get(s[1]);
//				System.out.println(s[0] +"  " + entry1);
//				System.out.println(s[1] +"  " + entry2);
//				double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
			}
		}
		br3.close();
		
		BufferedReader br4 = new BufferedReader(new FileReader(new File(disgenefile)));
		String line4 = null;
		Set<String> mappedset = new HashSet<>();
		
		while ((line4 = br4.readLine()) != null) {
			String [] s = line4.split("\\|");
//			System.out.println(s[0]);
			String disid = s[0];
			String gene = s[1];
			if (entry_index.containsKey(disid) && entry_index.containsKey(gene)) {
				mappedset.add(disid);
				int entry1 = entry_index.get(disid);
				int entry2 = entry_index.get(gene);
//					System.out.println(s[0] +"  " + entry1);
//					System.out.println(s[1] +"  " + entry2);
//					double score = Integer.parseInt(s[2]);
				cgg.addEdge(entry1, entry2);
			}
			
			// Additional mapping according to biological meaning
			
//			if (!entry_index.containsKey(disid) && extra) {
//				for (String id : dcn_omim.keySet()) {
//					List<String> omim_ids = dcn_omim.get(id);
//					if (omim_ids.contains(disid)) {
//						disid = id; // update id
//						if (entry_index.containsKey(disid) && entry_index.containsKey(gene)) {
//							int entry1 = entry_index.get(disid);
//							int entry2 = entry_index.get(gene);
//							cgg.addEdge(entry1, entry2);
//						}
//					}
//			}
			
			if(extra) {
				for (String id : dcn_omim.keySet()) {
					List<String> omim_ids = dcn_omim.get(id);
					if (!entry_index.containsKey(disid) && omim_ids.contains(disid)) {
						disid = id; // update id
						if (entry_index.containsKey(disid) && entry_index.containsKey(gene)) {
							int entry1 = entry_index.get(disid);
							int entry2 = entry_index.get(gene);
							cgg.addEdge(entry1, entry2);
						}
					}
				}
			}
		}
		br4.close();
		
		System.out.println("Total mapped disease: " + mappedset.size());
		
		for (int i=0; i<n; i++) {
			cgg.addEdge(i, i);
		}
		return cgg;
	}
	
	public static DisGraph createCommGeneticsGraph(DisGraph dg, DisGraph pp, 
			Map<Integer, List<Integer>> dis_genelist) {
		int N1 = dg.getNodes();
		int N2 = pp.getNodes();
		DisGraph cgg = new DisGraph(N1 + N2);
		
		//add disease comorbidity network
		for (int i=0; i<N1; i++) {
			Set<Integer> neibors = dg.getNeibor(i).getKeys();
			for (int neibor : neibors) {
				cgg.addEdge(i, neibor);
			}
		}
		
		//add protein-protein interaction network
		for (int i=0; i<N2; i++) {
			Set<Integer> neibors = pp.getNeibor(i).getKeys();
			for (int neibor : neibors) {
				int idx1 = i + N1;
				int idx2 = neibor + N1;
				cgg.addEdge(idx1, idx2);
			}
		}
		
		//add disease-gene connection
		for (int dis_idx: dis_genelist.keySet()) {
			for (int gene_idx: dis_genelist.get(dis_idx)) {
				cgg.addEdge(dis_idx, gene_idx + N1);
			}
		}
		
		
		return cgg;
	}
	
	
	/**
	 * Creates a UMLS-name mapping
	 * @param mapfiles mapping file
	 * @throws IOException
	 */
	public static void createIdNameMap(List<String> mapfiles) throws IOException {
		
		for (String mapfile : mapfiles) {
			BufferedReader br = new BufferedReader(new FileReader(mapfile));
			String line = null;
			while((line = br.readLine()) != null){
				String[] parts = line.split("\\|");
				idnamemap.put(parts[0], parts[1].toLowerCase());
				nameidmap.put(parts[1].toLowerCase(), parts[0]);
			}
			br.close();
		}
		
	}
	
	/**
	 * Create UMLS-name mapping for OMIM
	 * @param mapfiles
	 * @throws IOException
	 */
	public static void createOMIMIdNameMap(List<String> mapfiles) throws IOException {
		
		for (String mapfile : mapfiles) {
			BufferedReader br = new BufferedReader(new FileReader(mapfile));
			String line = null;
			while((line = br.readLine()) != null){
				String[] parts = line.split("\\|");
				OMIMidnamemap.put(parts[1], parts[0].toLowerCase());
				OMIMnameidmap.put(parts[0].toLowerCase(), parts[1]);
			}
			br.close();
		}
		
	}
	
	/**
	 * Write the bipartite graph to a file, format is "node1|node2|score"
	 * @param cgg the bipartite file
	 * @param filename a file to be written
	 * @throws IOException
	 */
	private static void writeCommGeneticsNet(DisGraph cgg, String filename) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
		for (int i = 0; i < cgg.getNet().size(); i++) {
			Set<Integer> sv = cgg.getNet().getKey(i);
			for (int k : sv) {
				String d1 = index_entry.get(i);
				String d2 = index_entry.get(k);
				double score = cgg.getNet().get(i, k);
				bw.write(d1 + "|" + d2 + "|" + score + "|" + "\n");
			}
			
		}
		bw.close();
	}
	
	
	
	public static void main(String[] args) throws IOException {
		
		
// 		#########################################################################
//		----------------------- Method 1: Build network from files  --------------------
//		#########################################################################
		
//		/*
//		 * Create additional DCN_OMIM map for give disease
//		 */
//		
//		String DCN_dis = "dementia";
//		String OMIM_pat = "alzheimer";
//		
//		DCNOMIMUMLSIDmap domim = new DCNOMIMUMLSIDmap(DCN_dis, OMIM_pat);
//		Map<String, List<String>> dcn_omim = domim.getDCNOMIMUMLSIDmap();
//		
//		/*
//		 * Build the bipartite network
//		 */
////		String commnetfile = "/Users/zhengc/Projects/AD_comorbidity/data/fares_comm_net_conf_ISMB_final_public.txt";
////		String dmnfile = "/Users/zhengc/Projects/AD_comorbidity/data/dmn_dm.txt";
////		String ppifile = "/Users/zhengc/Projects/AD_comorbidity/data/gene_gene_string_cut.txt";
////		String disgenefile = "/Users/zhengc/workspace/FARES_final/data/OMIM/mapped_OMIM/OMIM_disease_gene_umls_id_diso";
//		
//		String commnetfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_lift_final_abbr.txt";;
////		String dmnfile = "/Users/zhengc/workspace/FARES_final/analysis/AD_comorbidity/data/dmn_dm.txt";
//		String ppifile = "/Users/zhengc/workspace/FARES_final/analysis/AD_comorbidity/data/gene_gene_string_cut.txt";
//		String disgenefile = "/Users/zhengc/workspace/FARES_final/data/OMIM/mapped_OMIM/OMIM_disease_gene_umls_id_diso";
//		
//		
//		DisGraph cgg = null;
//		try {
////			cgg = CommGeneticsLoader.createCommGeneticsGraph(commnetfile, dmnfile, ppifile, disgenefile, dcn_omim);
//			cgg = CommGeneticsLoader.createCommGeneticsGraph(commnetfile, ppifile, disgenefile, dcn_omim);
//			System.out.println("Total edges: " + cgg.getEdges());
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		
		
		
		/* Write commgenetics file */
//		String netfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_PPI/DCN_PPI_net_lift.txt";
//		writeCommGeneticsNet(cgg, netfile);
//		
//		String disUMLS_name_file = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_PPI/disUMLS_name_lift.txt";
//		writeDCNIdName(commnetfile, disUMLS_name_file);
		
//		String dcnmapfile = "/Users/zhengc/workspace/FARES_final/data/FARES/map/term_umls_id_diso";
//		util.DCNOMIMUMLSIDmap.createDCNIdNameMap(dcnmapfile);
//		String diseaseid = util.DCNOMIMUMLSIDmap.dcnnameidmap.get(DCN_dis);
//		System.out.println(diseaseid);
//		int disease_idx = entry_index.get(diseaseid);
//		System.out.println(disease_idx);
//		System.out.println(cgg.getNeibor(disease_idx));
//		
//		Map<String, Double> commDisMap = ComSearch.searchFaers2(diseaseid, cgg);
//		System.out.println("\nNeibor node method......");
//		System.out.println("Number of comorbidities found: " + commDisMap.size() + "\n");
//		for (String dis : commDisMap.keySet()) {
//			System.out.println(dis);
//		}
		
		
// 		#########################################################################
//		--------------------- Method 2: Build network from  network -------------
//		#########################################################################
		
		String dcnmapfile = "/Users/zhengc/workspace/FARES_final/data/FARES/map/term_umls_id_diso";
		String omimmapfile = "/Users/zhengc/workspace/FARES_final/data/OMIM/map/OMIM_umls_id_diso";
		util.DCNOMIMUMLSIDmap.createDCNIdNameMap(dcnmapfile);
		util.DCNOMIMUMLSIDmap.createOMIMIdNameMap(omimmapfile);
		
		String DCN_dis = "dementia";
		String OMIM_pat = "alzheimer";
		
		DCNOMIMUMLSIDmap domim = new DCNOMIMUMLSIDmap(DCN_dis, OMIM_pat);
		Map<String, List<String>> dcn_omim = domim.getDCNOMIMUMLSIDmap();
		for (String dis: dcn_omim.keySet()) {
			System.out.println(dcn_omim.get(dis));
		}
		
		/*
		 * Build heterogenous network
		 */

		String rulefile = "/Users/zhengc/workspace/FARES_final/data/FARES/UMLS_map_data/"
				+ "ID_indications_all_clean_width_umls_id_diso_filtered_sp_lift_1_0.000002_3_processed.txt";
		String commnetfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_lift_final_abbr.txt";;
//		String dmnfile = "/Users/zhengc/workspace/FARES_final/analysis/AD_comorbidity/data/dmn_dm.txt";
		String ppifile = "/Users/zhengc/workspace/FARES_final/analysis/AD_comorbidity/data/gene_gene_string_cut.txt";
		String disgenefile = "/Users/zhengc/workspace/FARES_final/data/OMIM/mapped_OMIM/OMIM_disease_gene_umls_id_diso";
		
		CommLoader.readAssocRules(rulefile);
		DisGraph commnet = CommLoader.createCommNet(rulefile);
		CommLoader.readPPI(ppifile);
		DisGraph ppinet = CommLoader.createPPINet(ppifile);
		
		Map<Integer, List<Integer>> dis_gene = CommLoader.readDisGene(disgenefile, dcn_omim);
		
		
		DisGraph cgg2 = createCommGeneticsGraph(commnet, ppinet, dis_gene);
		System.out.println(cgg2.getEdges());
		System.out.println(cgg2.getNodes());
	}
}
