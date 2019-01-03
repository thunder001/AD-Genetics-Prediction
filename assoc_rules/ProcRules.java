package assoc_rules;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This class provides methods for processing mined association rules
 * @author zhengc
 *
 */
public class ProcRules {

	public static void main(String[] args) throws IOException {

		/* Command line arguments*/
		String rules = args[0];
		String out = args[1];
		String category = args[2];
		
		BufferedReader br = new BufferedReader(new FileReader(new File(rules)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out)));
		
		// Writer header line
		String header = category + "_1" + "|" + "Support_1 " + "|" + category + "_2" + "|" +"Support_2" + "|" +
		"Confidence" + "|" + "Lift" + "|" + "Leverage" + "|" + "Conviction" + "\n";
		bw.write(header);
		
		/* The following codes parse association rules from weka FPGrowth output 
		 * into easy handled format for downstream analysis 
		 * It takes several steps
		 */
		String line = br.readLine();
		//System.out.println(line);
		line = br.readLine();
		//System.out.println(line);
		while ((line = br.readLine()) != null) {
			String procline = "";
			
			//Step 1: Remove meaningless symbols including "=T", "<>", and count number in the begging etc
			//String pat = "=T|\\<?[a-z]+:|\\>|:|=|\\(|\\)|\\s*\\d+\\. +";
			String pat = "=1|\\<?[a-z]+:|\\>|:|=|\\(|\\)|\\s*\\d+\\. +";
			Pattern p = Pattern.compile(pat);
			Matcher m = p.matcher(line);
			line = m.replaceAll("");
			
			//Step 2: Distinguish spaces inside features from outside by replace spaces inside with "_"   
//			String pat0 = "(\\[([A-Z\\-\\d]+ +)+[A-Z\\-\\d]+\\])|(\\[([A-Z \\-\\d]+, )+([A-Z \\-\\d])+\\])";
			String pat0 = "(\\[\\s+\\w+\\]|\\[([\\w\\-\\s/]+\\s+)+[\\w\\-/]+\\])|(\\[([\\w\\s\\-/]+,\\s)+([\\w\\s\\-/])+\\])";
			Pattern p0 = Pattern.compile(pat0);
			Matcher m0 = p0.matcher(line);
			StringBuffer sb = new StringBuffer();
//			System.out.println(m0.find());
			while (m0.find()) {
				String replace = line.substring(m0.start(), m0.end()).replace(" ", "_");
//				System.out.println(replace);
				m0.appendReplacement(sb, replace);
			}
			m0.appendTail(sb);
			line = sb.toString();
			
			//Step 3: Replace spaces outside of features with "|"
			String delim = " +";
			String[] parts = line.split(delim);
			for (int i=0; i < parts.length; i++) {
				procline = String.join("|", procline, parts[i]);
			}
			procline = procline.replaceFirst("\\|", "");
			
			//Step 4: Change back spaces inside of features
			String pat4 = "_";
			Pattern p4 = Pattern.compile(pat4);
			Matcher m4 = p4.matcher(procline);
			procline = m4.replaceAll(" ");
			System.out.println(procline);
			
			bw.write(procline + "\n"); // Write to file
					
		}
		br.close();
		bw.close();

	}

}
