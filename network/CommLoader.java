package network;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
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

import GraphAlgorithm.GraphGenerator;
import GraphAlgorithm.RandGraph;

/**
 * Provides methods for building a disese comorbidity network (DCN) from FAERS \
 * and save it into a file
 * 
 * @author zhengc
 *
 */
public class CommLoader {
	public static HashSet<String> disease_set = new HashSet<String>();
	public static HashMap<String, Integer> disease_index = new HashMap<String, Integer>();
	public static HashMap<Integer, String> index_disease = new HashMap<Integer, String>();
	public static Map<String, String> idnamemap = new HashMap<String, String>();
	public static Map<String, String> nameidmap = new HashMap<String, String>();
	public static HashSet<String> soc_set = new HashSet<String>();
	public static Map<String, String> umlssocmap = new HashMap<String, String>();
	public static Map<String, String> umlssocabbrmap = new HashMap<String, String>();
	public static Map<String, Integer> soc_idx = new HashMap<String, Integer>();
	public static Map<Integer, String> idx_soc = new HashMap<Integer,String>();
	public static final double P = 0.85;
	
	public static HashSet<String> gene_set = new HashSet<String>();
	public static List<String> gene_list = new ArrayList<String>();
	public static HashMap<String, Integer> gene_index = new HashMap<String, Integer>();
	public static HashMap<Integer, String> index_gene = new HashMap<Integer, String>();
	
	
	
	/**
	 * Reads an association rule file to create necessary map file for building DCN
	 * @param filename an association rule file 
	 * @throws IOException
	 */
	public static void readAssocRules(String filename) throws IOException{
		// get all unique disease set
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			line = line.replace("[", "");
			line = line.replace("]", "");
//			System.out.println(line);
			String[] symptoms = line.split("\\|")[0].split(",");
//			System.out.println(symptoms.length);
			String[] diseases = line.split("\\|")[2].split(",");
			
			for (int i = 0; i < symptoms.length; i++) {
				symptoms[i] = symptoms[i].trim();
				disease_set.add(symptoms[i]);
			}
			for (int j = 0; j < diseases.length; j++) {
				diseases[j] = diseases[j].trim();
				disease_set.add(diseases[j]);
			}
		}
		br.close();
		
		int m = 0; 
		for (String d: disease_set) {
				disease_index.put(d, m);
				index_disease.put(m, d);
				m++;
		}
	}
	
	public static void readUMLSSOC(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		String line = null;
		while ((line=br.readLine()) != null) {
			String[] parts = line.split("\\$");
			String umls = parts[0];
			String soc = parts[8];
			String soc_abbr = parts[9];
			soc_set.add(soc);
			umlssocmap.put(umls, soc);
			umlssocabbrmap.put(umls, soc_abbr);
		}
		br.close();
		
		int m = 0;
		for (String soc : soc_set) {
			soc_idx.put(soc, m);
			idx_soc.put(m, soc);
			m++;
		}
		
	}
	
	
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
		
		int m = 0; 
		for (String d: gene_set) {
				gene_index.put(d, m);
				index_gene.put(m, d);
				m++;
		}
		
	}
	
	
	public static DisGraph createPPINet(String ppifile) throws IOException {
		File file = new File(ppifile);
		BufferedReader br = new BufferedReader(new FileReader(file));
		
		int n = gene_index.size();
		DisGraph ppinet = new DisGraph(n);
		
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			String[] parts = line.split("\\|");
			String dis1 = parts[0];
			int d1 = gene_index.get(dis1);
//			System.out.println(symptoms.length);
			String dis2 = parts[1];
			int d2 = gene_index.get(dis2);
			ppinet.addEdge(d1, d2);
		}
		br.close();
		return ppinet;
	}
	
	
	public static Map<Integer, List<Integer>> readDisGene(String disgenefile, 
			Map<String, List<String>> dcn_omim) throws IOException{
		
		Map<Integer, List<Integer>> dis_gene = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(new File(disgenefile)));
		String line = null;
		
		while ((line = br.readLine()) != null) {
			String [] s = line.split("\\|");
//			System.out.println(s[0]);
			String dis = s[0];
			String gene = s[1];
			if (disease_index.containsKey(dis)) {
				int dis_idx = disease_index.get(dis);
				if (gene_index.containsKey(gene)) {
					int gene_idx = gene_index.get(gene);
					if (!dis_gene.containsKey(dis_idx)) {
						List<Integer> genes = new ArrayList<Integer>();
						genes.add(gene_idx);
						dis_gene.put(dis_idx, genes);
					} else {
						dis_gene.get(dis_idx).add(gene_idx);
					}
				}
			} else {
				for (String id : dcn_omim.keySet()) {
					List<String> omim_ids = dcn_omim.get(id);
					if (omim_ids.contains(dis)) {
						dis = id; // update id
//						System.out.println(dis);
						if (disease_index.containsKey(dis) && gene_index.containsKey(gene)) {
							int dis_idx = disease_index.get(dis);
							int gene_idx = gene_index.get(gene);
							if (!dis_gene.containsKey(dis_idx)) {
								List<Integer> genes = new ArrayList<Integer>();
								genes.add(gene_idx);
								dis_gene.put(dis_idx, genes);
							} else {
								dis_gene.get(dis_idx).add(gene_idx);
							}
						}
					}
				}
			}
			
			// add additional dcn_omim map
			
			
			
//			for (String id : dcn_omim.keySet()) {
//				List<String> omim_ids = dcn_omim.get(id);
//				if (!disease_index.containsKey(dis) && omim_ids.contains(dis)) {
//					dis = id; // update id
//					if (disease_index.containsKey(dis) && gene_index.containsKey(gene)) {
//						int dis_idx = disease_index.get(dis);
//						int gene_idx = gene_index.get(gene);
//						dis_gene.get(dis_idx).add(gene_idx);
//					}
//				}
//			}
		}
		br.close();
		return dis_gene;
	} 
	
	/**
	 * Builds a DCN from an association rule file
	 * @param filename an association rule file
	 * @return an undirected and unweighted DisGraph
	 * @throws IOException
	 */
	public static DisGraph createCommNet(String filename) throws IOException{
		
		File file = new File(filename);
		BufferedReader br = new BufferedReader(new FileReader(file));
		
		int n = disease_index.size();
		DisGraph commnet = new DisGraph(n);
		
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
		
			line = line.replace("[", "");
			line = line.replace("]", "");
		
			//undirected, unweighted network 
		 
			String[] disease1 = line.split("\\|")[0].split(",");
			List<String> diseases_raw = new ArrayList<String>(Arrays.asList(disease1));
			String[] disease2 = line.split("\\|")[2].split(",");
			diseases_raw.addAll(Arrays.asList(disease2));
			
			List<String> diseases = new ArrayList<String>();
			for (String d:diseases_raw) {
				d = d.trim();
				diseases.add(d);
			}
			
			for (int i=0; i<diseases.size() - 1; i++) {
				//System.out.println(diseases.get(i));
				int d1 = disease_index.get(diseases.get(i));
				for (int j=i+1; j<diseases.size(); j++) {
					int d2 = disease_index.get(diseases.get(j));
					//double conf = Double.parseDouble(line.split("\\|")[4]);
					//commnet.addNode(d1, d2, conf);
					commnet.addEdge(d1, d2);
				}
			}
			
		}
		br.close();
		
		return commnet;
	}
	
	
	public static DisGraph createRandNet(String filename) throws IOException {
		DisGraph commnet = createCommNet(filename);
		DisGraph randnet = (new RandGraph(commnet)).getRandGraph();
		return randnet;
	}
	
	
	public static DisGraph createERNet(String filename) throws IOException {
		DisGraph commnet = createCommNet(filename);
		int V = commnet.getNodes();
		double p = (double) V / (V*(V-1)/2);
		System.out.println("Probability: " + p);
		DisGraph randnet = GraphGenerator.simple(V, p);
		return randnet;
	}
	
	public static DisGraph createCommNetFromAdjFile(String netfile, String adjfile) throws IOException {
		DisGraph commnet = createCommNet(netfile);
		
		DisGraph adjnet = new DisGraph(commnet.getNodes());
		
		BufferedReader br = new BufferedReader(new FileReader(new File(adjfile)));
		String line = null;
		int num_line = 0;
		while ((line=br.readLine()) != null) {
			String[] parts = line.split("\\|");
			for (int i=0; i<parts.length; i++) {
				int entry = Integer.parseInt(parts[i]);
				if (entry==1) {
					adjnet.addEdge(num_line, i);
				}
			}
			num_line++;
		}
		br.close();
		return adjnet;
	}
	
	
	public static List<String> getNeibors(String dis, DisGraph dg) {
		
		String dis_id = CommLoader.nameidmap.get(dis);
		int dis_idx = CommLoader.disease_index.get(dis_id);
		
		SparseVector sv = dg.getNeibor(dis_idx);
		List<String> neibors = new ArrayList<>();
		for (int j : sv.getKeys()) {
			String id = CommLoader.index_disease.get(j);
			neibors.add(CommLoader.idnamemap.get(id));
		}
		Collections.sort(neibors);
		return neibors;
	}
	
	public static Map <Integer, Integer> getDegDistri(DisGraph dg) {
		Map <Integer, Integer> dis_neibor = new HashMap<>();
		SparseMatrix commMatrix = dg.getNet();
		for (int i=0; i<commMatrix.size(); i++) {
			dis_neibor.put(i, commMatrix.getNeibor(i).nnz());
		}
		
		return dis_neibor;
	}
	
	public static List<String> compareGraph(DisGraph dg, DisGraph randgraph,  String disname1, String disname2) {
		
		List<String> neibors_original = getNeibors(disname1, dg);
		
		System.out.println("\nComorbidities of " + disname1 + " in original graph:" + neibors_original.size() + "\n");
		for (String d : neibors_original) {
			System.out.println(d);
		}
		
		
		List<String> neibors_random = getNeibors(disname2, randgraph);
		System.out.println("\nComorbidities of " + disname2 + " in random graph:" + neibors_random.size() + "\n" );
		for (String d : neibors_random) {
			System.out.println(d);
		}
		
		List<String> neibors_common = new ArrayList<>();
		if (neibors_original.size() >= neibors_random.size()) {
			for (String j: neibors_original) {
				if (neibors_random.contains(j)) {
					neibors_common.add(j);
				}
			}
		} else {
			for (String j: neibors_random) {
				if (neibors_original.contains(j)) {
					neibors_common.add(j);
				}
			}
		}
		
		Collections.sort(neibors_common);
		System.out.println("\nShared comorbidities of " + disname1 + "and " + disname2 + 
				" between original and random graphs :" + neibors_common.size() + "\n" );
		for (String d : neibors_common) {
			System.out.println(d);
		}
		return neibors_common;
	}
	
	
	/**
	 * Write network to a txt file
	 * @param omimnet a DisGraph
	 * @param filename a file to be written
	 * @throws IOException
	 * dis1_UMLS| dis1_name|dis1_SOC|dis1_SOC_idx|dis2_UMLS|dis2_name|dis2_SOC|dis2_SOC_idx|conf
	 */
	
	public static void writeCommNet(DisGraph commnet, String filename) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
		bw.write("#dis1_UMLS" + "|" + "dis1_name" + "|"+ "dis1_SOC" + "|" + "dis1_SOC_abbr" + "|" +  "dis1_SOC_idx" + "|" + 
		"dis2_UMLS" + "|" + "dis2_name" + "|"+ "dis2_SOC" + "|" +  "dis2_SOC_abbr" + "|" + "dis2_SOC_idx" + "|" + "conf" + "\n");
		for (int i = 0; i < commnet.getNet().size(); i++) {
			String d1 = index_disease.get(i);
			String d1_name = idNameConverter(d1, idnamemap);
			String d1_soc = umlssocmap.get(d1);
			String d1_soc_abbr = umlssocabbrmap.get(d1);
			int d1_soc_idx = 0;
			if (soc_idx.keySet().contains(d1_soc)) {
				d1_soc_idx = soc_idx.get(d1_soc);
			} else {
				d1_soc_idx = soc_idx.size();
			}
			
			Set<Integer> sv = commnet.getNet().getKey(i);
			for (int k : sv) {
				String d2 = index_disease.get(k);
				String d2_name = idNameConverter(d2, idnamemap);
				String d2_soc = umlssocmap.get(d2);
				String d2_soc_abbr = umlssocabbrmap.get(d2);
				int d2_soc_idx = 0;
				if (soc_idx.keySet().contains(d2_soc)) {
					d2_soc_idx = soc_idx.get(d2_soc);
				} else {
					d2_soc_idx = soc_idx.size();
				}
				
				double conf = commnet.getNet().get(i, k) ;
				bw.write(d1 + "|" + d1_name + "|"+ d1_soc + "|" + d1_soc_abbr + "|" + d1_soc_idx + "|" 
				+ d2 + "|" + d2_name + "|"+ d2_soc + "|" + d2_soc_abbr + "|" + d2_soc_idx + "|" + conf + "\n");
			}
		}
		bw.close();
	}
	

	
	/**
	 * Write disease node (UMLS) to a txt file
	 * @param filename a txt file to be written
	 * @throws IOException
	 */
	private static void writeTerms(String filename) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
		for (String dis : disease_index.keySet()) {
			bw.write(dis + "\n");
		}
		bw.close();
	}
	
	/**
	 * Write degree for all disease nodes to a file
	 * @param filename a file to be written
	 * @throws IOException
	 */
	public static void writeNodeDeg(DisGraph commnet, String filename) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
		Map <Integer, Integer> dis_neibor = getDegDistri(commnet);
		
		bw.write("disease" + "|" + "numofdeg" + "index" + "\n");;
		for (int k : dis_neibor.keySet()) {
			String dis = idNameConverter(index_disease.get(k), idnamemap);
			int deg = dis_neibor.get(k);
			bw.write(dis + "|" + deg + "|" + k + "\n");
		}
		bw.close();
	}

	/**
	 * Converts a disease UMLS ID to a disease concept name
	 * @param id a UMLS ID
	 * @param idnamemap a Map from UMLS to name
	 * @return a disease name
	 */
	private static String idNameConverter(String id, Map<String, String> idnamemap) {
		return idnamemap.get(id);
	}

	/**
	 * Creates a map from UMLS to disease concept name
	 * @param mapfile a file containing map information from UMLS to disease concept name
	 * @throws IOException
	 */
	public static void createIdNameMap(String mapfile) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(mapfile));
		String line = null;
		while((line = br.readLine()) != null){
			String[] parts = line.split("\\|");
			idnamemap.put(parts[0], parts[1].toLowerCase());
			nameidmap.put(parts[1].toLowerCase(), parts[0]);
		}
		br.close();
	}
	
			
	/*
	 * main method	
	 */
	
	public static void main(String[] args) throws IOException {
		/* Create disease comorbidity network
		 * Write network to a file
		 * Query network
		 */
//		**************** rule file *************************
		
		String rulefile = "/Users/zhengc/workspace/FARES_final/data/FARES/UMLS_map_data/"
				+ "ID_indications_all_clean_width_umls_id_diso_filtered_sp_lift_1_0.000002_3_processed.txt";

//		**************** UMLS id to disease name map file file *************************

		String mapfile = "/Users/zhengc/workspace/FARES_final/data/FARES/map/umls_id_name_diso";		
		
//		**************** network file to be saved *************************
		String netfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_conf_final_abbr.txt";
//		String netfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_lift_final_abbr.txt";

//		**************** term file to be saved *************************
		String termfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_term_lift_final_public.txt";	

		
		
//		####################################################################
//		------------------ Build and save network to files ----------------
//		####################################################################
		
		readAssocRules(rulefile);
		readUMLSSOC("/Users/zhengc/workspace/FARES_final/analysis/Cluster/umls_soc.txt");

		DisGraph commnet = createCommNet(rulefile);

		//		System.out.println(commnet.getNet());
		System.out.println("Nodes: " + commnet.getNodes());
		System.out.println("Edges: " + commnet.getEdges());
		

		createIdNameMap(mapfile);
 		
//		writeCommNet(commnet, netfile);
//		writeTerms(termfile);
		
		
		// Print all diseases
//		List<String> disList = new ArrayList<String>();
//		for (String dis : disease_index.keySet()) {
//			dis.trim();
//			disList.add(idnamemap.get(dis));
//		}
//		Collections.sort(disList);
//		
//		System.out.println();
//		System.out.println(disList.size());
//		for (String dis : disList) {
//			System.out.println(dis);
//		}

		
		// display neibor degree distribution
//		Map <Integer, Integer> dis_neibor = getDegDistri(commnet);
//		
//		
//		System.out.println("Number of disease comorbidity distribution: ");
//		for (int k : dis_neibor.keySet()) {
//			System.out.println(k + ":" + idnamemap.get(index_disease.get(k)) + ": " + dis_neibor.get(k));
//		}
//		
//		String degreefile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_degree_lift_adj.txt";
//		writeNodeDeg(commnet, degreefile);
		

		
		
//		for (int j=0; j<200; j++) {
//			System.out.println(j + ":\n" + commMatrix.getNeibor(j));
//		}
		
		
//		####################################################################
//		---------------------Randomize network analysis -----------------------------
//		####################################################################
		
//		String adjfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_random_deg/net_adj_2.txt";
//		String netdir = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_random_er/";
//
//		// String netdir = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_random_er/";
////		DisGraph rdg = createRandNet(rulefile);
//		DisGraph rdg = createERNet(rulefile);
////		DisGraph rdg = createCommNetFromAdjFile(rulefile, adjfile);
//		
//		// compare random graph with original graph: comorbidities for specific disease
//		String dis1 = "dementia";
//		String dis2 = "dementia";
//		List<String> neibors_comm = compareGraph(commnet, rdg, dis1, dis2);
//		
//		String rnetfile = netdir + "lift_5" + ".txt";
//		CommLoader.writeCommNet(rdg,rnetfile);

		

//		####################################################################
//		------------------ Build and save randomized networks to files ----------------
//		####################################################################		
		
//		String netdir = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_random/";
//		int N = 1000;
//		for (int i=0; i<N; i++) {
//			System.out.printf("\nRandom graph %d\n",i);
//			
//			DisGraph rg = createRandNet(rulefile);
//			System.out.println(rg.getEdges());
//			
//			String rnetfile = netdir + "lift_" + i + ".txt";
//			CommLoader.writeCommNet(rg,rnetfile);
//			
//			String rdegreefile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN_random/fares_comm_net_degree_lift.txt";
//			CommLoader.writeNodeDeg(rg, rdegreefile);
//
//		}
		
		
	}
}	


