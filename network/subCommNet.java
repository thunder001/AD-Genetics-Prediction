package network;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class subCommNet {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		
//		String disfile = "/Users/zhengc/Projects/AD_comorbidity/data/AD_commobidity.txt";
//		String netfile = "/Users/zhengc/workspace/FARES/data/FARES/network/fares_comm_net_conf_ISMB_final_public.txt";
//		String subnetfile = "/Users/zhengc/Projects/AD_comorbidity/data/fares_comm_net_sub_AD.txt";
		
		String disfile = "/Users/zhengc/workspace/FARES_final/analysis/AD_comorbidity/data/AD_commobidity_rw_lift.txt";
		String netfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_lift_final_abbr.txt";
		String subnetfile = "/Users/zhengc/workspace/FARES_final/analysis/network/DCN/fares_comm_net_lift_sub_AD.txt";
		List<String> dislist = getDisList(disfile);
		for (String dis: dislist) {
			System.out.println(dis);
		}
		createSubCommNet(dislist, netfile, subnetfile);
	}
	
//	public CommGraph createSubCommNet(List<String> dlist, CommGraph cg) {
//		SparseMatrix comMatrix = cg.getNet();
//		CommGraph subcommnet = new CommGraph(dlist.size());
//		for (String dis : dlist) {
//			SparseVector neibors = cg.getNet().getNeibor(dis);
//			
//		}
//		
//		return subcommnet;
//	}
	public static List<String> getDisList(String disfile) throws IOException{
		List<String> dislist = new ArrayList<>();
		BufferedReader br = new BufferedReader(new FileReader(new File(disfile)));
		String line = null;
		while((line = br.readLine()) != null) {
			String dis = line.split("\\|")[0];
			dislist.add(dis);
		}
		br.close();
		return dislist;
	}
	
	public static void createSubCommNet(List<String> dlist, String netfile, String subnetfile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(netfile)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(subnetfile)));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			String[] parts = line.split("\\|");
			String dis1 = parts[1];
			String dis2 = parts[6].replace("-- ", "");
			if (dlist.contains(dis1) & dlist.contains(dis2)) {
				line = line.replace("-- ", "");
				bw.write(line + "\n");
			}
		}
		br.close();
		bw.close();
		
	}

}
