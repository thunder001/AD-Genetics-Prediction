package GraphAlgorithm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import network.CommGeneticsLoader;
import network.CommLoader;
import network.DisGraph;
import network.SparseVector;
import util.DCNOMIMUMLSIDmap;

public class DeNovoPred {
	
	private Map<String, Integer> generank = new HashMap<>();
	
	/**
	 * Initializes a new DovovoPrediction object
	 * @param cgg a CommGeneticsGraph
	 * @param dis the disease whose genes will be de novo predicted 	
	 * @throws IOException if dis doesn't exist
	 */
	public DeNovoPred(DisGraph cgg, String dis, double p) throws IOException{
		// check if the disease is in disease comorbidity network
		String dis_id = null;
		if (DCNOMIMUMLSIDmap.dcnnameidmap.containsKey(dis)) {
			dis_id = DCNOMIMUMLSIDmap.dcnnameidmap.get(dis);
		} else {
			System.out.println("No such disease in this network!");
		}
		// Critical part for removing all disease gene edges
		 int dis_idx = CommGeneticsLoader.entry_index.get(dis_id); // for CommGeneticsloader from netfile
		// int dis_idx = CommLoader.disease_index.get(dis_id); // for CommGeneticsloader from net
		System.out.printf("Disease: %s\n", dis);
		System.out.printf("Disease index: %s\n", dis_idx);
		System.out.printf("Disease id: %s\n", dis_id);
		System.out.println("Disease_gene edges before removing: " + cgg.getNet().getNeibor(dis_idx));
		cgg.removeEdges(dis_idx); // key code for do novo prediction
		System.out.println("Disease_gene edges after removing: " + cgg.getNet().getNeibor(dis_idx));
		
		// Random walk to rank disease genes
		List<String> dislist = new ArrayList<>();
		dislist.add(dis_id);
		SparseVector seed = RandomWalk.createSeedVector(dislist);
		RandomWalk rm = new RandomWalk(cgg, seed, p);
		generank = rm.getRWRank();
	}
	
	/**
	 * Get de novo prediction gene ranks for specific disease
	 * @return de novo prediction gene ranks for specific disease
	 */
	public Map<String, Integer> getDoNovo() {
		return generank;
	}
	
	/**
	 * Write gene-rank pairs into a CSV file, format is "Gene,Rank,Percentage"
	 * @param dn de novo prediction map
	 * @param dnfile a file to be written
	 * @throws IOException
	 */
	public void saveDoNovo(Map<String, Integer> dn, String dnfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dnfile)));
		bw.write("Gene" + ","  + "Rank" + "," + "Percentage" +  "\n");
		for (String gene : dn.keySet()) {
//			System.out.println(sortedgenemap.get(s));
			int rank = dn.get(gene);
			double perc = rank * 1.0 / CommGeneticsLoader.num_gene * 100;
			bw.write(gene + "," + rank + "," + perc + "\n");
		}
		bw.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		// step 1: build DCN-PPI network
		System.out.println("Building DCN network...");
		
//		String commnetfile = args[0];
//		String ppifile = args[1];
//		String disgenefile = args[2];
		
		String commnetfile = "./data/fares_comm_net_lift_final_abbr.txt";
		String ppifile = "./data/gene_gene_string_cut.txt";
		String disgenefile = "./data/OMIM_disease_gene_umls_id_diso";
		
		
//		String dcnmapfile = args[3];
//		String omimmapfile = args[4];
		
		String dcnmapfile = "./data/term_umls_id_diso";
		String omimmapfile = "./data/OMIM_umls_id_diso";
		util.DCNOMIMUMLSIDmap.createDCNIdNameMap(dcnmapfile);
		util.DCNOMIMUMLSIDmap.createOMIMIdNameMap(omimmapfile);
		String DCN_dis = "dementia";
		String OMIM_pat = "alzheimer";
		DCNOMIMUMLSIDmap domim = new DCNOMIMUMLSIDmap(DCN_dis, OMIM_pat);
		Map<String, List<String>> dcn_omim = domim.getDCNOMIMUMLSIDmap();
		
		DisGraph cgg = null;
		try {
			cgg = CommGeneticsLoader.createCommGeneticsGraph(commnetfile, ppifile, disgenefile, dcn_omim);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("Done!");
		
		
		// step 2: optimize restart probability
//		String rulefile = args[5];
//		String donovoDir = args[6]
//		String predType = args[7]
				
		String rulefile = "./data/ID_indications_all_clean_width_umls_id_diso_filtered_sp_lift_1_0.000002_3_processed.txt";
		String donovoDir = "./results/evaluation/";
		String predType = "pred" //"optimize", "random" or "pred"
		
		if (preType=="optimize") {
			/**
			 * do novo prediction for different restart probabilities
			 */
			System.out.println("Starting do novo gene prediction for " + DCN_dis + "!");
			List<Double> ps = new ArrayList<Double>();
			for(int i=1; i<10; i++){
				ps.add(i * 0.1);
			}
			
			int count = 1;
			for(double p : ps){
				System.out.println("\nTesting " + count + " of " + ps.size());
				DeNovoPred dnp = new DeNovoPred(cgg, DCN_dis, p);
				Map<String, Integer> generank = dnp.getDoNovo();
				String donovofile = donovoDir + "commnet_opti_" + P + DCN_dis  + ".csv";
				dnp.saveDoNovo(generank, donovofile);
				count++;
			}
		}
		
		// step 3: get gene ranks for randomized network
		if (preType=="random") {
//	 		#########################################################################
//			----------------------- Method: Erdos-Renyi graph ---------------------------------
//			#########################################################################
			System.out.println("\nStaring random node do novo prediction......\n");
		
//			/*
//			 * Build the bipartite network
	
			CommLoader.readAssocRules(rulefile);
			DisGraph commnet = CommLoader.createCommNet(rulefile);
			int V = commnet.getNodes();
			double p = (double) V / (V*(V-1)/2);
			System.out.println("Probability: " + p);
			
			CommLoader.readPPI(ppifile);
			DisGraph ppinet = CommLoader.createPPINet(ppifile);
			
			List<String> entry_list = new ArrayList<String>();
			List<String> disease_list = new ArrayList<String>();
			disease_list = util.CollectionsEx.setToList(CommLoader.disease_set);
			entry_list.addAll(disease_list);
			entry_list.addAll(CommLoader.gene_list);
			CommGeneticsLoader.createIndex(entry_list);
			CommGeneticsLoader.num_disease = disease_list.size();
			CommGeneticsLoader.entry_list = entry_list;
			CommGeneticsLoader.num_gene = CommLoader.gene_list.size();
			
			
			Map<Integer, List<Integer>> dis_gene = CommLoader.readDisGene(disgenefile, dcn_omim);
			
			int N = 1000;
			double P = 0.5;
			
			for (int i=0; i<N; i++) {
				System.out.println();
				if (i > 10 && i % 100 == 0) {
					System.out.printf("Finish %d of %d prediction!", i, N);
				}
				
				DisGraph randnet = GraphGenerator.simple(V, p);
				System.out.println("Total nodes in randamized disease graph: " + randnet.getNodes());
				System.out.println("Total edges in randomized disease graph: " + randnet.getEdges());
				DisGraph cgg = CommGeneticsLoader.createCommGeneticsGraph(randnet, ppinet, dis_gene); 
				System.out.println("Total nodes in randamized heterogeneous graph: " + cgg.getNodes());
				System.out.println("Total edges in randomized heterogeneous graph: " + cgg.getEdges());
				DeNovoPred dnp = new DeNovoPred(cgg, DCN_dis, P); //!speed limit operation
				Map<String, Integer> generank = dnp.getDoNovo();
				
				System.out.println("Size of gene rank:" + generank.size());
				
				
				String donovofile = donovoDir + "random_er_" + P + "_" + i + "_" + DCN_dis  + ".csv";
				dnp.saveDoNovo(generank, donovofile);

			}
		}
		
		//Step 4: get gene rank from DCN
		if (preType=="pred") {
			DeNovoPred dnp = new DeNovoPred(cgg, DCN_dis, 0.5); 
			Map<String, Integer> generank = dnp.getDoNovo();
			
			System.out.println("Size of gene rank:" + generank.size());
			String donovofile = donovoDir + "commnet_" + 0.5 + "_" + DCN_dis  + ".csv";
			dnp.saveDoNovo(generank, donovofile);
		}
		
	}

}
