package network;

import java.util.Iterator;
import java.util.Map;

/**
 * Provides routine methods for DisGraph operation
 * @author zhengc
 *
 */
public class DisGraph {

	private SparseMatrix disnet;

	/**
	 * Initialize a DisGraph represented as a SparseMatrix
	 * @param n number of nodes
	 */
	public DisGraph(int n) {
		disnet = new SparseMatrix(n);
	}

	/**
	 * Add edge between two nodes (Integer)
	 * @param d1
	 * @param d2
	 * @param score
	 */
	public void addEdge(int d1, int d2, double score) {
		
		disnet.put(d1, d2, score);
		disnet.put(d2, d1, score);
	}
	
	/**
	 * Add edge between two nodes (Integer), weight is set to 1
	 * @param d1
	 * @param d2
	 */
	public void addEdge(int d1, int d2) {
		
		disnet.put(d1, d2, 1);
		disnet.put(d2, d1, 1);
	}
	
	/**
	 * Add edge between two nodes (CommNode)
	 * @param node1
	 * @param node2
	 */
	public void addEdge(CommNode node1, CommNode node2) {
		
		int d1 = node1.getIndex();
		int d2 = node2.getIndex();
		disnet.put(d1, d2, 1);
		disnet.put(d2, d1, 1);
	}

	/**
	 * Get neighbor nodes for given node 
	 * @param d
	 * @return a SparseVector representing the neighbor nodes
	 */
	
	public SparseVector getNeibor(int d) {
		return disnet.getNeibor(d);
	}
	
	
	
	/**
	 * Remove one edge
	 * @param d1 node 1
	 * @param d2 node 2
	 */
	public void removeEdge(int d1, int d2) {
		
		disnet.put(d1, d2, 0);
		disnet.put(d2, d1, 0);
	}
	
	/**
	 * Remove all edges between disease and its OMIM genes
	 * @param dis a disease node
	 */
	public void removeEdges(int dis) {
		
		SparseVector nb = disnet.getNeibor(dis);
		
		Iterator<Integer> keys = nb.iterator();
		while(keys.hasNext()) {
			int key = keys.next();
//			System.out.println(key);
			if (key > CommGeneticsLoader.num_disease) {
//				System.out.println("Removing link " + key);
				keys.remove();
				removeEdge(dis, key);
			}
		}
	}
	

	/**
	 * Get number of edges of the DisGraph
	 * @return number of edges of the DisGraph
	 */
	public int getEdges() {
		return this.getNet().nnz() / 2; // symmetric network
	}

	/**
	 * Get number of nodes of the DisGraph
	 * @return number of nodes of the DisGraph
	 */
	public int getNodes() {
		return disnet.size();
	}

	/**
	 * Get the SparseMatrix of the DisGraph
	 * @return
	 */
	public SparseMatrix getNet() {
	
			return disnet;
		}

}