package assoc_rules;

import java.io.File;
import java.io.FileWriter;
import java.time.LocalDateTime;

import weka.associations.Apriori;
import weka.associations.FPGrowth;
import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;

/**
 * This class implements FP growth algorithm to obtain association rule from clean, UMLS mapped and filtered FAERS data
 * @author zhengc
 *
 */
public class AssocFPGrowth{
	public static void main(String args[]) throws Exception{
		//load dataset
		
		String dataset = args[0];
		System.out.println("Algorithm starts: " + LocalDateTime.now());
		System.out.println();
		
		DataSource source = new DataSource(dataset);
		//get instances object 
		Instances data = source.getDataSet();
		int numAttr = data.numAttributes();
		int numIns =  data.numInstances();
		
		System.out.println("Number of attributes is " + numAttr);
		System.out.println("Number of instances is " + numIns);
		
		// lift-based algorithm
		String[] opts = new String[10];
		opts[0] = "-T"; opts[1] = "1";
		opts[2] = "-C"; opts[3] = "1";
		opts[4] = "-M"; opts[5] = "0.000002";
		opts[6] = "-N"; opts[7] ="200000";
		opts[8] = "-I"; opts[9] = "3";
		
		model.setOptions(opts);
		//build model
		model.buildAssociations(data);
		
		//print out and save the extracted rules to file
		String results = model.toString();
		System.out.println(results);

		
		FileWriter fw = new FileWriter(new File(resultfile));
		fw.write(results);
		fw.close();
		
		System.out.println("Algorithm ends: " + LocalDateTime.now());
	}
}
