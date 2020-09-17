package dataprep;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.Normalizer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class LongToWide {
	
	public static void main(String[] args) throws IOException {
		
		Map<Integer, Set<String>> pid_cids = new HashMap<>();
		
		String lfname1 = "/Users/zhengc/Projects/AD_metaformin/data/medical_history/processed/ADNI1_umls_id_dsyn_2.csv";
//		String lfname2 = "/Users/zhengc/Projects/AD_metaformin/data/medical_history/processed/ADNI3_umls_id_dsyn_2.csv";
		String lfname3 = "/Users/zhengc/Projects/AD_metaformin/data/drug/processed/drug_umls_id.csv";
		
		Map<Integer, Set<String>> pid_cids1 = getCids(lfname1);
//		Map<Integer, Set<String>> pid_cids2 = getCids(lfname2);
		Map<Integer, Set<String>> pid_cids3 = getCids(lfname3);
		
		Set<Integer> pids = pid_cids1.keySet();
		Set<Integer> pid3 = pid_cids3.keySet();
		pids.retainAll(pid3);
		
		for (int pid: pids) {
			Set<String> cids = pid_cids1.get(pid);
			Set<String> cids3 = pid_cids3.get(pid);
			cids.addAll(cids3);
			pid_cids.put(pid, cids);
		}
		
		
		
		String wfname = "/Users/zhengc/Projects/AD_metaformin/data/drug/processed/Comb_drug_mh_umls_id_width.csv";
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(wfname)));
		
		for (int pid: pid_cids.keySet()) {
			Set<String> cids = pid_cids.get(pid);
			StringBuilder sb = new StringBuilder();
			for (String cid: cids) {
				sb.append(cid + ",");
			}
			String cidstr0 = sb.toString();
			String cidstr = cidstr0.replaceAll("[,]$", "");
//			String cidstr = cidstr0.substring(0, cidstr0.length() - 1);
			bw.write(cidstr + "\n");
		}
		
		bw.close();
		
	}
	
	private static Map<Integer, Set<String>> getCids(String lfname) throws IOException {
		Map<Integer, Set<String>> pid_cids = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(new File(lfname)));
		String line = br.readLine();
		while ((line=br.readLine()) != null) {
			String[] parts = line.split(",");
			System.out.println(parts.length);
			int pid = Integer.parseInt(parts[0]);
			System.out.println(pid);
			String cid = parts[2];
			System.out.println(cid);
			String ad_status = null;
			if (parts.length > 3) {
				ad_status = parts[3];
			}
			
			if (!pid_cids.containsKey(pid)) {
				Set<String> cids = new HashSet<>();
				cids.add(cid);
				if (ad_status != null && (ad_status.equals("AD") || ad_status.equals("MCI"))) cids.add("C0002395");
				pid_cids.put(pid, cids);
				
			} else {
				pid_cids.get(pid).add(cid);
				if (ad_status != null && (ad_status.equals("AD") || ad_status.equals("MCI"))) pid_cids.get(pid).add("C0002395");
			}
		}
		br.close();
		return pid_cids;
	}
}
	
	
	

