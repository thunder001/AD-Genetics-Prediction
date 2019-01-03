package GraphAlgorithm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import network.CommGeneticsLoader;
import network.DisGraph;
import network.SparseMatrix;
import network.SparseVector;
import util.DCNOMIMUMLSIDmap;
import util.FileToList;

/**
 * Implement random walk with restart algorithm using sparse matrix data structure
 * @author zhengc
 *
 */
public class RandomWalk {
	
//	static final double P = 0.15;
//	static final double P = 0.15;
	static final double err = 0.000001;
	private Map<String, Integer> generank;
	
	/**
	 * Constructs a new RandomWalk object from a graph and seed vector
	 * @param cgg a CommGeneticsGraph 
	 * @param E the seed vector
	 * @param p restart rate
	 */
	public RandomWalk(DisGraph cgg, SparseVector E, double p) {
		Map<String, Double> genescore = new HashMap<>();
		int size = cgg.getNodes();
		SparseVector x = new SparseVector(size);
		SparseMatrix A = cgg.getNet();
		A = normalizeMatrix(A);
		int iter = 0;
		SparseVector last_x = E;
		while((last_x.plus(x.scale(-1)).norm()) > err) {
			last_x = x;
			x = A.times(x).scale(1-p);
			x = x.plus(E.scale(p));
//			System.out.print(x);
			iter++;
//			System.out.println("Finish iteration " + iter + ":");
		}
		System.out.println("Total iteration " + iter);
		
		genescore = svToMap(x, E);
		generank = scoreToRank(genescore);
		
	}
	
	/**
	 * Get the gene-rank pairs from random walk with restart
	 * @return the gene-rank pairs from random walk with restart
	 */
	public  Map<String, Integer> getRWRank() {
		return generank;
	}
	
	
	/**
	 * Generate column/left stochastic matrix by normalization
	 * @param comMatrix
	 * @return normalized SparseMatrix
	 */
	private static SparseMatrix normalizeMatrix(SparseMatrix comMatrix) {
		
		int n = comMatrix.size();
		double[] total_w = new double[n]; //Total weight of a column
		SparseMatrix A = new SparseMatrix(n);
		
		for (int j = 0; j < n; j++){
			total_w[j] = 0;
			for (int i = 0; i < n; i++){
				total_w[j] += comMatrix.get(i, j);
			}
		}

		
		//normalize the adjacency matrix
		for (int i = 0; i < n; i++){
			
			for (int j = 0; j < n; j++){
				double norm_w = comMatrix.get(i, j);
				norm_w = norm_w * 1.0 / total_w[j];
				A.put(i, j, norm_w);
			}
		}
		
		return A;
	}
	
	/**
	 * Create a seed from a list
	 * @param nodesids a list of seed nodes
	 * @return the seed vector
	 * @throws IOException
	 */
	public static SparseVector createSeedVector(List<String> nodesids) throws IOException{
		
		int n = CommGeneticsLoader.entry_index.size();
		SparseVector E = new SparseVector(n);

		for (String id : nodesids) {
			if (CommGeneticsLoader.entry_index.containsKey(id)) {
				int index = CommGeneticsLoader.entry_index.get(id);
//				E.put(index, 1);
				E.put(index, 1.0 / nodesids.size());
			} else {
				System.out.println(id + " is NOT included in our network!");
			}
		}
		return E;
	}
	
	
	public static SparseVector createSeedVectorFromIdx(List<Integer> nodeidxs) {
		int n = CommGeneticsLoader.entry_index.size();
		SparseVector E = new SparseVector(n);

		for (int idx : nodeidxs) {
			if (idx > n - 1) {
				System.out.println(idx + " is out of network!");
			} else {
				E.put(idx, 1.0 / nodeidxs.size());;
			}
		}
		return E;
	}
	
	/**
	 * Convert a SparseVector from random walk into a HashMap reversely ordered by value
	 * @param sv the SparseVector representing gene rank score result
	 * @param seeds_vector the seed vector
	 * @return a HashMap reversely ordered by value 
	 */
	private static Map<String, Double> svToMap(SparseVector sv, SparseVector seeds_vector) {
		Map<String, Double> genemap = new HashMap<String,Double>();
//		System.out.println(CommGeneticsLoader.num_disease);
//		System.out.println(CommGeneticsLoader.entry_list.size());
		for (int i = CommGeneticsLoader.num_disease; i < CommGeneticsLoader.entry_list.size(); i++){
//			double val_seed = seeds_vector.get(i);
			Set<Integer> seeds_key = seeds_vector.getKeys();
			if (!seeds_key.contains(i)) { // Test if the node is one of seeds
				double val = sv.get(i);
//				System.out.println(index_entry.get(num_disease));
				
				if (val != 0){
					String gene = CommGeneticsLoader.index_entry.get(i);
					genemap.put(gene, val);
				}
			}
		}
		genemap = util.CollectionsEx.sortByValueR(genemap);
		return genemap;
	}
	
	/**
	 * Convert a gene rank score into a gene rank HashMap reversely ordered by value
	 * @param sortedgenemap a reversely sorted gene score map
	 * @return a gene rank  HashMap
	 */
	private static Map<String, Integer> scoreToRank(Map<String, Double> sortedgenemap) {
		Map<String, Integer> generank = new HashMap<String,Integer>();
		int rank = 0;
		for (String s : sortedgenemap.keySet()) {
//			System.out.println(sortedgenemap.get(s));
			rank++;
			generank.put(s, rank);
			generank = util.CollectionsEx.sortByValue(generank);
		}
		
		return generank;
	}
	
//	System.out.println(genemap.size());
	
	
	
	/**
	 * Write gene rank scores into CSV file, file format is "Rank, Gene, Score"
	 * @param sortedgenemap HashMap sorted by descend values
	 * @param rankfile file to be written
	 * @throws IOException
	 */
	public static void saveScorefile(Map<String, Double> sortedgenemap, String rankfile) throws IOException {
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(rankfile)));
		int rank = 0;
		bw.write("Rank" + "," + "Gene" + "," + "Score" + "," + "\n");
		for (String s : sortedgenemap.keySet()) {
//			System.out.println(sortedgenemap.get(s));
			Double score = sortedgenemap.get(s);
			rank++;
			bw.write(rank + "," + s + "," + score + "," + "\n");
		}
		bw.close();
	}
	
	/**
	 * Write gene rank file into a CSV file, file format is "Rank, Gene"
	 * @param generank gene rank HashMap
	 * @param rankfile the file to be written
	 * @throws IOException
	 */
	public static void saveRankfile(Map<String, Integer> generank, String rankfile) throws IOException {
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(rankfile)));
		bw.write("Rank" + "," + "Gene" + ","  + "\n");
		for (String s : generank.keySet()) {
//			System.out.println(sortedgenemap.get(s));
			int rank = generank.get(s);
			bw.write(rank + "," + s + "," + "\n");
		}
		bw.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		
		/*
		 * Create additional DCN_OMIM map for give disease
		 */
		String dcnmapfile = "./data/term_umls_id_diso";
		String omimmapfile = "./data/OMIM_umls_id_diso";
		util.DCNOMIMUMLSIDmap.createDCNIdNameMap(dcnmapfile);
		util.DCNOMIMUMLSIDmap.createOMIMIdNameMap(omimmapfile);
		
		String DCN_dis = "dementia";
		String OMIM_pat = "alzheimer";
		
		DCNOMIMUMLSIDmap domim = new DCNOMIMUMLSIDmap(DCN_dis, OMIM_pat);
		Map<String, List<String>> dcn_omim = domim.getDCNOMIMUMLSIDmap();
		
		/*
		 * Build the bipartite network
		 */

		String commnetfile = "./results/fares_comm_net_lift_final_abbr.txt";
		String ppifile = "./data/gene_gene_string_cut.txt";
		String disgenefile = "./data/OMIM_disease_gene_umls_id_diso";
		DisGraph cgg = null;
		try {
			cgg = CommGeneticsLoader.createCommGeneticsGraph(commnetfile, ppifile, disgenefile, dcn_omim);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/*
		 * Create seed vector
		 */
		int n = cgg.getNodes();
		
		List<String> seeds = new ArrayList<>();
		String seed_dis1 = "dementia";
		String seed_dis1_umls = null;
		if (util.DCNOMIMUMLSIDmap.dcnnameidmap.containsKey(seed_dis1)) {
			seed_dis1_umls = util.DCNOMIMUMLSIDmap.dcnnameidmap.get(seed_dis1);
			System.out.println(seed_dis1 + ": " + seed_dis1_umls);
			seeds.add(seed_dis1_umls);
		} else {
			System.out.println(seed_dis1 + " is NOT included in our network!");
		}
		
		String seed_genes_file = "./data/AD_omim_genes.txt";
		List<String> seed_genes = FileToList.createListFromfile(seed_genes_file);
		
		seeds.addAll(seed_genes);
//		diseases.add(disease2);
		SparseVector seeds_vec = new SparseVector(n);
		try {
			seeds_vec = createSeedVector(seeds);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Seed vector: " + seeds_vec);
		
		/*
		 * Random walk
		 */

		RandomWalk rwrank = new RandomWalk(cgg, seeds_vec, 0.5);
		Map<String, Integer> result = rwrank.getRWRank();

		/*
		 * Process result and write to file
		 */

		String rankfile = "./results/AD_novel_genes.csv";
		try {
			saveRankfile(result, rankfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
