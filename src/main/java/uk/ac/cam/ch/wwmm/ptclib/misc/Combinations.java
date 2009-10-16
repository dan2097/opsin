package uk.ac.cam.ch.wwmm.ptclib.misc;

import java.util.ArrayList;
import java.util.List;
/** For a list of integers, generate all lists of non-negative integers where all
 * of the list members are strictly lower than their corresponding integer in the
 * input list.
 * 
 * @author ptc24
 *
 */
public final class Combinations {

	private List<List<Integer>> combList;
	private List<Integer> scheme;
	
	private Combinations() {
		
	}
	
	/**For a list of integers, generate all lists of non-negative integers where
	 * all of the list members are strictly lower than their corresponding
	 * integer in the input list.
	 * 
	 * @param scheme The input list of integers.
	 * @return The list of lists.
	 */
	public static List<List<Integer>> makeCombinations(List<Integer> scheme) {
		Combinations c = new Combinations();
		c.combList = new ArrayList<List<Integer>>();
		c.scheme = scheme;
		c.grow(new ArrayList<Integer>());
		return c.combList;
	}
	
	private void grow(List<Integer> soFar) {
		if(soFar.size() == scheme.size()) {
			combList.add(soFar);
		} else {
			for(int i=0;i<scheme.get(soFar.size());i++) {
				List<Integer> ll = new ArrayList<Integer>();
				ll.addAll(soFar);
				ll.add(i);
				grow(ll);
			}
		}
	}
	
}
