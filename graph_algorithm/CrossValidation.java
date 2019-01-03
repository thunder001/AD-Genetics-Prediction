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
import network.DisGraph;
import network.SparseVector;
import util.DCNOMIMUMLSIDmap;

/**
 * Implements leave one out cross validation to test the performance of disease genetics prediction
 * @author zhengc
 *
 */
public class CrossValidation {

	private Map<String, Integer> generank = new HashMap<>();
	
	/**
	 * Initiate and create a cross validation object 
	 * @param cgg CommGeneticsGraph
	 * @param dis query disease
	 * @param genelist omim genes for query disease
	 * @throws IOException
	 */
	public CrossValidation(DisGraph cgg, String dis, List<String> genelist) throws IOException {
		String dis_id = null;
		if (DCNOMIMUMLSIDmap.dcnnameidmap.containsKey(dis)) {
			dis_id = DCNOMIMUMLSIDmap.dcnnameidmap.get(dis);
		} else {
			System.out.println("No such disease in this network!");
		}
		
		List<String> dislist = new ArrayList<>();
		dislist.add(dis);
		SparseVector seed = RandomWalk.createSeedVector(dislist);
		int dis_idx = CommGeneticsLoader.entry_index.get(dis_id);
		
		int count = 1;
		for (String gene : genelist) {
			System.out.println(count + " of " +  genelist.size());
			if (CommGeneticsLoader.entry_index.containsKey(gene)) {
				int gene_idx = CommGeneticsLoader.entry_index.get(gene);
				cgg.removeEdge(dis_idx, gene_idx);
				RandomWalk rw = new RandomWalk(cgg, seed, 0.2);
				Map<String, Integer> result = rw.getRWRank();
				int rank = result.get(gene);
				generank.put(gene, rank);
				count++;
			}
		}
	}
	
	/**
	 * Get disease OMIM gene ranking from leave one out cross validation 
	 * @return disease OMIM gene ranking map
	 */
	public Map<String, Integer> getCV() {
		return generank;
	}
	
	/** 
	 * Write leave one out cross validation to a file
	 * @param cv disease OMIM gene ranking map
	 * @param cvfile a file to be written
	 * @throws IOException
	 */
	public static void saveCV(Map<String, Integer> cv, String cvfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(cvfile)));
		bw.write("Gene" + ","  + "Rank" + "," + "Percentage" +  "\n");
		for (String gene : cv.keySet()) {
//			System.out.println(sortedgenemap.get(s));
			int rank = cv.get(gene);
			double perc = rank / CommGeneticsLoader.num_gene * 100;
			bw.write(gene + "," + rank + "," + perc + "\n");
		}
		bw.close();
	}
	
	public static void main(String[] args) throws IOException {
		
		System.out.println("Building network...");
		
		String dcnmapfile = "/Users/zhengc/workspace/FARES/data/FARES/map/term_umls_id_diso";
		String omimmapfile = "/Users/zhengc/workspace/FARES/data/OMIM/map/OMIM_umls_id_diso";
		util.DCNOMIMUMLSIDmap.createDCNIdNameMap(dcnmapfile);
		util.DCNOMIMUMLSIDmap.createOMIMIdNameMap(omimmapfile);
		
		String DCN_dis = "dementia";
		String OMIM_pat = "alzheimer";
		
		DCNOMIMUMLSIDmap domim = new DCNOMIMUMLSIDmap(DCN_dis, OMIM_pat);
		Map<String, List<String>> dcn_omim = domim.getDCNOMIMUMLSIDmap();
		
		/*
		 * Build the bipartite network
		 */
		String commnetfile = "/Users/zhengc/Projects/AD_comorbidity/data/fares_comm_net_conf_ISMB_final_public.txt";
		String dmnfile = "/Users/zhengc/Projects/AD_comorbidity/data/dmn_dm.txt";
		String ppifile = "/Users/zhengc/Projects/AD_comorbidity/data/gene_gene_string_cut.txt";
		String disgenefile = "/Users/zhengc/workspace/FARES/data/OMIM/mapped_OMIM/OMIM_disease_gene_umls_id_diso";
		DisGraph cgg = null;
		try {
//			cgg = CommGeneticsLoader.createCommGeneticsGraph(commnetfile, dmnfile, ppifile, disgenefile, dcn_omim);
			cgg = CommGeneticsLoader.createCommGeneticsGraph(commnetfile, ppifile, disgenefile, dcn_omim);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("Done!");
		
		/*
		 * Cross validation
		 */
		System.out.println("Starting crossvalidate...");
		
		String genelistfile = "/Users/zhengc/Projects/AD_comorbidity/data/AD_omim_genes.txt";
		List<String> genelist = util.FileToList.createListFromfile(genelistfile);
		
		CrossValidation cv = new CrossValidation(cgg, DCN_dis, genelist);
		Map<String, Integer> generank = cv.getCV();
		
		for (String gene : generank.keySet()) {
			System.out.println(gene + ": " + generank.get(gene));
		}
		
		String cv_result_file = "/Users/zhengc/Projects/AD_comorbidity/results/evaluation/cv_result.csv";
		saveCV(generank, cv_result_file);

	}
}
