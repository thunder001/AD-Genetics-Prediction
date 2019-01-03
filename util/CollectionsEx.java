package util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Provides additional methods for collection operation
 * @author zhengc
 *
 */
public class CollectionsEx {
	@SuppressWarnings("hiding")
	/**
	 * Returns sorted HashMap by descend Value
	 * @param map the original HashMap 
	 * @return sorted HashMap by descend Value
	 */
	public static <String, Double extends Comparable<? super Double>> Map<String, Double> sortByValueR( Map<String, Double> map ) {
		/* sort map by value reversed*/
		List<Map.Entry<String, Double>> list = new LinkedList<>(map.entrySet());
		Collections.sort(list, Collections.reverseOrder(new Comparator<Map.Entry<String, Double>>() {
			
			public int compare( Map.Entry<String, Double> o1, Map.Entry<String, Double> o2 ) {
				return ( o1.getValue() ).compareTo( o2.getValue() );
			}
		}) );
		

		Map<String, Double> result = new LinkedHashMap<>();
		for (Map.Entry<String, Double> entry : list){
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
	}
	
	@SuppressWarnings("hiding")
	/**
	 * Returns sorted HashMap by Value
	 * @param map the original HashMap 
	 * @return sorted HashMap by Value
	 */
	public static <String, Double extends Comparable<? super Double>> Map<String, Double> sortByValue( Map<String, Double> map ) {
		/* sort map by value*/
		List<Map.Entry<String, Double>> list = new LinkedList<>(map.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String, Double>>() {
			
			public int compare( Map.Entry<String, Double> o1, Map.Entry<String, Double> o2 ) {
				return ( o1.getValue() ).compareTo( o2.getValue() );
			}
		} );
		

		Map<String, Double> result = new LinkedHashMap<>();
		for (Map.Entry<String, Double> entry : list){
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
	}
	
//	public static Map<String, String> sortByKey(Map<String, String> map ) {
//		/* sort map by key*/
//		Map<String, String> sortedmap = new HashMap<String, String>();
//		Set<String> keys = map.keySet();
//		List<String> key_list = setToList(keys);
//		Collections.sort(key_list);
//		
//		for(String key : key_list) {
//			System.out.println(key);
//			sortedmap.put(key, map.get(key));
//		}
//		return sortedmap;
//	}
	
	/**
	 * Creates a list from a set
	 * @param set a set
	 * @return a list
	 */
	public static List<String> setToList(Set<String> set) {
		List<String> list = new ArrayList<String>();
		for (String s: set) {
			list.add(s);
		}
		Collections.sort(list);
		return list;
	}

}
