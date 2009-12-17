package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.Map;

/**
 * Provides valency checking features and a lookup on the possible valencies
 * for an atom given its element and charge
 *
 * Also used to perform a final check on the output of OPSIN, to reject interpretations
 * that result in hypervalent structures due to incorrect names or misinterpreted names
 *
 * @author ptc24/dl387
 *
 */
class ValencyChecker {

	private static final Map<String, Integer> expectedDefaultValency;//used to decide on the likely valency state
	private static final Map<String, Integer> valencyInHW;//used to decide whether an atom has spare valency in a ring, these are the same as specified in the Hantzch-Widman system
	private static final Map<String, HashMap<Integer, Integer[]>> possibleStableValencies;//used to decide on the likely valency state

	static {
		expectedDefaultValency = new HashMap<String, Integer>();
		expectedDefaultValency.put("B", 3);
		expectedDefaultValency.put("Al", 3);
		expectedDefaultValency.put("In", 3);
		expectedDefaultValency.put("Ga", 3);
		expectedDefaultValency.put("Tl", 3);
		expectedDefaultValency.put("C", 4);
		expectedDefaultValency.put("Si", 4);
		expectedDefaultValency.put("Ge", 4);
		expectedDefaultValency.put("Sn", 4);
		expectedDefaultValency.put("Pb", 4);
		expectedDefaultValency.put("N", 3);
		expectedDefaultValency.put("P", 3);
		expectedDefaultValency.put("As", 3);
		expectedDefaultValency.put("Sb", 3);
		expectedDefaultValency.put("Bi", 3);
		expectedDefaultValency.put("O", 2);
		expectedDefaultValency.put("S", 2);
		expectedDefaultValency.put("Se", 2);
		expectedDefaultValency.put("Te", 2);
		expectedDefaultValency.put("Po", 2);
		expectedDefaultValency.put("F", 1);
		expectedDefaultValency.put("Cl", 1);
		expectedDefaultValency.put("Br", 1);
		expectedDefaultValency.put("I", 1);
		expectedDefaultValency.put("At", 1);

		//in order of priority in the HW system
		valencyInHW = new HashMap<String, Integer>();
		valencyInHW.put("F", 3);//IUPAC says 1 for halogens, but this makes no sense as in a ring valency will never be 1
		valencyInHW.put("Cl", 3);
		valencyInHW.put("Br", 3);
		valencyInHW.put("I", 3);
		valencyInHW.put("O", 2);
		valencyInHW.put("S", 2);
		valencyInHW.put("Se", 2);
		valencyInHW.put("Te", 2);
		valencyInHW.put("N", 3);
		valencyInHW.put("P", 3);
		valencyInHW.put("As", 3);
		valencyInHW.put("Sb", 3);
		valencyInHW.put("Bi", 3);
		valencyInHW.put("Si", 4);
		valencyInHW.put("Ge", 4);
		valencyInHW.put("Sn", 4);
		valencyInHW.put("Pb", 4);
		valencyInHW.put("B", 3);
		valencyInHW.put("Al", 3);
		valencyInHW.put("Ga", 3);
		valencyInHW.put("In", 3);
		valencyInHW.put("Tl", 3);
		valencyInHW.put("Hg", 2);

		possibleStableValencies = new HashMap<String, HashMap<Integer, Integer[]>>();
		possibleStableValencies.put("H", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("He", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Li", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Be", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("B", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("C", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("N", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("O", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("F", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ne", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Na", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Mg", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Al", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Si", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("P", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("S", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Cl", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ar", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("K", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ca", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ga", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ge", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("As", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Se", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Br", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Kr", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Rb", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Sr", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("In", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Sn", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Sb", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Te", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("I", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Xe", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Cs", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ba", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Tl", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Pb", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Bi", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Po", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("At", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Rn", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Fr", new HashMap<Integer, Integer[]>());
		possibleStableValencies.put("Ra", new HashMap<Integer, Integer[]>());

		possibleStableValencies.get("H").put(0, new Integer[]{1});
		possibleStableValencies.get("He").put(0, new Integer[]{0});
		possibleStableValencies.get("Li").put(0, new Integer[]{1});
		possibleStableValencies.get("Be").put(0, new Integer[]{2});
		possibleStableValencies.get("B").put(0, new Integer[]{3});
		possibleStableValencies.get("C").put(0, new Integer[]{4});
		possibleStableValencies.get("N").put(0, new Integer[]{3});
		possibleStableValencies.get("O").put(0, new Integer[]{2});
		possibleStableValencies.get("F").put(0, new Integer[]{1});
		possibleStableValencies.get("Ne").put(0, new Integer[]{0});
		possibleStableValencies.get("Na").put(0, new Integer[]{1});
		possibleStableValencies.get("Mg").put(0, new Integer[]{2});
		possibleStableValencies.get("Al").put(0, new Integer[]{3});
		possibleStableValencies.get("Si").put(0, new Integer[]{4});
		possibleStableValencies.get("P").put(0, new Integer[]{3,5});
		possibleStableValencies.get("S").put(0, new Integer[]{2,4,6});
		possibleStableValencies.get("Cl").put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get("Ar").put(0, new Integer[]{0});
		possibleStableValencies.get("K").put(0, new Integer[]{1});
		possibleStableValencies.get("Ca").put(0, new Integer[]{2});
		possibleStableValencies.get("Ga").put(0, new Integer[]{3});
		possibleStableValencies.get("Ge").put(0, new Integer[]{4});
		possibleStableValencies.get("As").put(0, new Integer[]{3,5});
		possibleStableValencies.get("Se").put(0, new Integer[]{2,4,6});
		possibleStableValencies.get("Br").put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get("Kr").put(0, new Integer[]{0,2});
		possibleStableValencies.get("Rb").put(0, new Integer[]{1});
		possibleStableValencies.get("Sr").put(0, new Integer[]{2});
		possibleStableValencies.get("In").put(0, new Integer[]{3});
		possibleStableValencies.get("Sn").put(0, new Integer[]{4});
		possibleStableValencies.get("Sb").put(0, new Integer[]{3,5});
		possibleStableValencies.get("Te").put(0, new Integer[]{2,4,6});
		possibleStableValencies.get("I").put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get("Xe").put(0, new Integer[]{0,2,4,6,8});
		possibleStableValencies.get("Cs").put(0, new Integer[]{1});
		possibleStableValencies.get("Ba").put(0, new Integer[]{2});
		possibleStableValencies.get("Tl").put(0, new Integer[]{3});
		possibleStableValencies.get("Pb").put(0, new Integer[]{4});
		possibleStableValencies.get("Bi").put(0, new Integer[]{3,5});
		possibleStableValencies.get("Po").put(0, new Integer[]{2,4,6});
		possibleStableValencies.get("At").put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get("Rn").put(0, new Integer[]{0,2,4,6,8});
		possibleStableValencies.get("Fr").put(0, new Integer[]{1});
		possibleStableValencies.get("Ra").put(0, new Integer[]{2});

		possibleStableValencies.get("H").put(1, new Integer[]{0});
		possibleStableValencies.get("Li").put(1, new Integer[]{0});
		possibleStableValencies.get("Be").put(2, new Integer[]{0});
		possibleStableValencies.get("B").put(1, new Integer[]{4});
		possibleStableValencies.get("B").put(-1, new Integer[]{2,4});
		possibleStableValencies.get("C").put(1, new Integer[]{3});
		possibleStableValencies.get("C").put(-1, new Integer[]{3});
		possibleStableValencies.get("N").put(1, new Integer[]{4});
		possibleStableValencies.get("N").put(-1, new Integer[]{2});
		possibleStableValencies.get("N").put(-2, new Integer[]{1});
		possibleStableValencies.get("O").put(1, new Integer[]{4});
		possibleStableValencies.get("O").put(1, new Integer[]{3});
		possibleStableValencies.get("O").put(-1, new Integer[]{1});
		possibleStableValencies.get("O").put(-2, new Integer[]{0});
		possibleStableValencies.get("F").put(1, new Integer[]{2});
		possibleStableValencies.get("F").put(-1, new Integer[]{0});
		possibleStableValencies.get("Na").put(1, new Integer[]{0});
		possibleStableValencies.get("Na").put(-1, new Integer[]{0});
		possibleStableValencies.get("Mg").put(2, new Integer[]{0});
		possibleStableValencies.get("Al").put(3, new Integer[]{0});
		possibleStableValencies.get("P").put(1, new Integer[]{4});
		possibleStableValencies.get("P").put(-1, new Integer[]{2});
		possibleStableValencies.get("S").put(1, new Integer[]{3});
		possibleStableValencies.get("S").put(-1, new Integer[]{1});
		possibleStableValencies.get("S").put(-2, new Integer[]{0});
		possibleStableValencies.get("Cl").put(1, new Integer[]{2});
		possibleStableValencies.get("Cl").put(-1, new Integer[]{0});
		possibleStableValencies.get("K").put(1, new Integer[]{0});
		possibleStableValencies.get("K").put(-1, new Integer[]{0});
		possibleStableValencies.get("Ca").put(2, new Integer[]{0});
		possibleStableValencies.get("Ga").put(3, new Integer[]{0});
		possibleStableValencies.get("Ga").put(-1, new Integer[]{2});
		possibleStableValencies.get("Ge").put(4, new Integer[]{0});
		possibleStableValencies.get("As").put(1, new Integer[]{4});
		possibleStableValencies.get("As").put(-1, new Integer[]{2});
		possibleStableValencies.get("As").put(-2, new Integer[]{1});
		possibleStableValencies.get("As").put(-3, new Integer[]{0});
		possibleStableValencies.get("Se").put(1, new Integer[]{3});
		possibleStableValencies.get("Se").put(-1, new Integer[]{1});
		possibleStableValencies.get("Se").put(-2, new Integer[]{0});
		possibleStableValencies.get("Br").put(1, new Integer[]{2});
		possibleStableValencies.get("Br").put(-1, new Integer[]{0});
		possibleStableValencies.get("Rb").put(1, new Integer[]{0});
		possibleStableValencies.get("Rb").put(-1, new Integer[]{0});
		possibleStableValencies.get("Sr").put(2, new Integer[]{0});
		possibleStableValencies.get("In").put(3, new Integer[]{0});
		possibleStableValencies.get("Sn").put(4, new Integer[]{0});
		possibleStableValencies.get("Sn").put(-1, new Integer[]{3});
		possibleStableValencies.get("Sb").put(3, new Integer[]{0});
		possibleStableValencies.get("Sb").put(1, new Integer[]{4});
		possibleStableValencies.get("Te").put(1, new Integer[]{3});
		possibleStableValencies.get("Te").put(-1, new Integer[]{1});
		possibleStableValencies.get("Te").put(-2, new Integer[]{0});
		possibleStableValencies.get("I").put(1, new Integer[]{2});
		possibleStableValencies.get("I").put(-1, new Integer[]{0});
		possibleStableValencies.get("Cs").put(1, new Integer[]{0});
		possibleStableValencies.get("Cs").put(-1, new Integer[]{0});
		possibleStableValencies.get("Ba").put(2, new Integer[]{0});
		possibleStableValencies.get("Pb").put(-1, new Integer[]{3});
		possibleStableValencies.get("Pb").put(2, new Integer[]{0});
		possibleStableValencies.get("Bi").put(3, new Integer[]{0});
		possibleStableValencies.get("Bi").put(1, new Integer[]{4});
		possibleStableValencies.get("At").put(-1, new Integer[]{0});
		possibleStableValencies.get("Fr").put(1, new Integer[]{0});
		possibleStableValencies.get("Ra").put(2, new Integer[]{0});
	}

	/**
	 * Given an element symbol (e.g. Na) and charge (e.g. 1) returns the highest stable valency that OPSIN knows is possible
	 * If for the particular combination of element symbol and charge the highest stable valency is not known null is returned
	 * @param symbol
	 * @param charge
	 * @return
	 */
	private static Integer getMaxValency(String symbol, int charge) {
		if (possibleStableValencies.get(symbol)!=null){
			if (possibleStableValencies.get(symbol).get(charge)!=null){
				return possibleStableValencies.get(symbol).get(charge)[possibleStableValencies.get(symbol).get(charge).length-1];
			}
		}
		return null;
	}

	/**
	 * Checks whether the total incoming valency to an atom exceeds its expected valency
	 * @param a
	 * @return
	 */
	static boolean checkValency(Atom a) {
		int valency = a.getIncomingValency();
		Integer maxVal;
		if (a.getValency()!=null){
			maxVal=a.getValency();
		}
		else{
			String symbol = a.getElement();
			int charge = a.getCharge();
			maxVal=getMaxValency(symbol, charge);
			if(maxVal==null) return true;
		}
		return valency <= maxVal;
	}

	/** Check whether valency is available on the atom to form a bond of the given order.
	 * spareValency and outValency are not taken into account.
	 * @param a atom you are interested in
	 * @param bondOrder order of bond required
	 */
	static boolean checkValencyAvailableForBond(Atom a, int bondOrder) {
		int valency =a.getIncomingValency() +bondOrder;
		Integer maxVal;
		if (a.getValency()!=null){
			maxVal=a.getValency();
		}
		else{
			String symbol = a.getElement();
			int charge = a.getCharge();
			maxVal = getMaxValency(symbol, charge);
			if(maxVal==null) return true;
		}
		return valency <= maxVal;
	}

	/** Check whether changing to a heteroatom will result in valency being exceeded
	 * spareValency and outValency is taken into account
	 * @param a atom you are interested in
	 * @param heteroatom element symbol of heteroatom
	 */
	static boolean checkValencyAvailableForReplacementByHeteroatom(Atom a, String heteroatom) {
		int valency =a.getIncomingValency();
		valency +=a.hasSpareValency() ? 1 : 0;
		valency +=a.getOutValency();
		Integer maxValOfHeteroAtom = getMaxValency(heteroatom, 0);
        return maxValOfHeteroAtom == null || valency <= maxValOfHeteroAtom;
		}

	/**
	 * Returns the default valency of an element or null if unknown
	 * @param element
	 * @param charge
	 * @return
	 */
	static Integer getDefaultValency(String element, int charge) {
		if (charge != 0){return null;}
		return expectedDefaultValency.get(element);
	}

	/**
	 * Returns the valency of an element in the HW system (useful for deciding whether something should have double bonds in a ring) or null if unknown
	 * Note that the HW system makes no claim about valency when the atom is charged
	 * @param element
	 * @return
	 */
	static Integer getHWValency(String element) {
		return valencyInHW.get(element);
	}

	/**
	 * Returns the maximum valency of an element or null if unknown
	 * @param element
	 * @param charge
	 * @return
	 */
	static Integer getMaximumValency(String element, int charge) {
		return getMaxValency(element, charge);
	}

	/**
	 * Returns the maximum valency of an element or null if unknown
	 * @param element
	 * @param charge
	 * @return
	 */
	static Integer[] getPossibleValencies(String element, int charge) {
		if (possibleStableValencies.get(element)==null){return null;}
		return possibleStableValencies.get(element).get(charge);
	}
}
