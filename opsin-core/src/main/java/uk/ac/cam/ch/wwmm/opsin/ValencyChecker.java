package uk.ac.cam.ch.wwmm.opsin;

import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;

/**
 * Provides valency checking features and a lookup on the possible valencies
 * for an atom given its element and charge
 *
 * Also used to perform a final check on the output of OPSIN, to reject interpretations
 * that result in hypervalent structures due to incorrect names or misinterpreted names
 *
 * @author ptc24
 * @author dl387
 *
 */
class ValencyChecker {

	/** used to decide on the likely valency state*/
	private static final Map<ChemEl, Integer> expectedDefaultValency = new EnumMap<ChemEl, Integer>(ChemEl.class);
	
	/** used to decide whether an atom has spare valency in a ring, these are the same as specified in the Hantzch-Widman system */
	private static final Map<ChemEl, Integer> valencyInHW = new EnumMap<ChemEl, Integer>(ChemEl.class);
	
	/** used to decide on the likely valency state */
	private static final Map<ChemEl, Map<Integer, Integer[]>> possibleStableValencies = new EnumMap<ChemEl, Map<Integer, Integer[]>>(ChemEl.class);

	static {
		expectedDefaultValency.put(ChemEl.B, 3);
		expectedDefaultValency.put(ChemEl.Al, 3);
		expectedDefaultValency.put(ChemEl.In, 3);
		expectedDefaultValency.put(ChemEl.Ga, 3);
		expectedDefaultValency.put(ChemEl.Tl, 3);
		expectedDefaultValency.put(ChemEl.C, 4);
		expectedDefaultValency.put(ChemEl.Si, 4);
		expectedDefaultValency.put(ChemEl.Ge, 4);
		expectedDefaultValency.put(ChemEl.Sn, 4);
		expectedDefaultValency.put(ChemEl.Pb, 4);
		expectedDefaultValency.put(ChemEl.N, 3);
		expectedDefaultValency.put(ChemEl.P, 3);
		expectedDefaultValency.put(ChemEl.As, 3);
		expectedDefaultValency.put(ChemEl.Sb, 3);
		expectedDefaultValency.put(ChemEl.Bi, 3);
		expectedDefaultValency.put(ChemEl.O, 2);
		expectedDefaultValency.put(ChemEl.S, 2);
		expectedDefaultValency.put(ChemEl.Se, 2);
		expectedDefaultValency.put(ChemEl.Te, 2);
		expectedDefaultValency.put(ChemEl.Po, 2);
		expectedDefaultValency.put(ChemEl.F, 1);
		expectedDefaultValency.put(ChemEl.Cl, 1);
		expectedDefaultValency.put(ChemEl.Br, 1);
		expectedDefaultValency.put(ChemEl.I, 1);
		expectedDefaultValency.put(ChemEl.At, 1);

		//in order of priority in the HW system
		valencyInHW.put(ChemEl.F, 1);
		valencyInHW.put(ChemEl.Cl, 1);
		valencyInHW.put(ChemEl.Br, 1);
		valencyInHW.put(ChemEl.I, 1);
		valencyInHW.put(ChemEl.O, 2);
		valencyInHW.put(ChemEl.S, 2);
		valencyInHW.put(ChemEl.Se, 2);
		valencyInHW.put(ChemEl.Te, 2);
		valencyInHW.put(ChemEl.N, 3);
		valencyInHW.put(ChemEl.P, 3);
		valencyInHW.put(ChemEl.As, 3);
		valencyInHW.put(ChemEl.Sb, 3);
		valencyInHW.put(ChemEl.Bi, 3);
		valencyInHW.put(ChemEl.Si, 4);
		valencyInHW.put(ChemEl.Ge, 4);
		valencyInHW.put(ChemEl.Sn, 4);
		valencyInHW.put(ChemEl.Pb, 4);
		valencyInHW.put(ChemEl.B, 3);
		valencyInHW.put(ChemEl.Al, 3);
		valencyInHW.put(ChemEl.Ga, 3);
		valencyInHW.put(ChemEl.In, 3);
		valencyInHW.put(ChemEl.Tl, 3);
		valencyInHW.put(ChemEl.Hg, 2);

		possibleStableValencies.put(ChemEl.H, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.He, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Li, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Be, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.B, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.C, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.N, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.O, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.F, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ne, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Na, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Mg, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Al, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Si, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.P, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.S, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Cl, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ar, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.K, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ca, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ga, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ge, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.As, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Se, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Br, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Kr, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Rb, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Sr, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.In, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Sn, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Sb, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Te, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.I, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Xe, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Cs, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ba, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Tl, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Pb, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Bi, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Po, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.At, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Rn, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Fr, new HashMap<Integer, Integer[]>());
		possibleStableValencies.put(ChemEl.Ra, new HashMap<Integer, Integer[]>());

		possibleStableValencies.get(ChemEl.H).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.He).put(0, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Li).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Be).put(0, new Integer[]{2});
		possibleStableValencies.get(ChemEl.B).put(0, new Integer[]{3});
		possibleStableValencies.get(ChemEl.C).put(0, new Integer[]{4});
		possibleStableValencies.get(ChemEl.N).put(0, new Integer[]{3});
		possibleStableValencies.get(ChemEl.O).put(0, new Integer[]{2});
		possibleStableValencies.get(ChemEl.F).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Ne).put(0, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Na).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Mg).put(0, new Integer[]{2});
		possibleStableValencies.get(ChemEl.Al).put(0, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Si).put(0, new Integer[]{4});
		possibleStableValencies.get(ChemEl.P).put(0, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.S).put(0, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Cl).put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Ar).put(0, new Integer[]{0});
		possibleStableValencies.get(ChemEl.K).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Ca).put(0, new Integer[]{2});
		possibleStableValencies.get(ChemEl.Ga).put(0, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Ge).put(0, new Integer[]{4});
		possibleStableValencies.get(ChemEl.As).put(0, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Se).put(0, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Br).put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Kr).put(0, new Integer[]{0,2});
		possibleStableValencies.get(ChemEl.Rb).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Sr).put(0, new Integer[]{2});
		possibleStableValencies.get(ChemEl.In).put(0, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Sn).put(0, new Integer[]{2,4});
		possibleStableValencies.get(ChemEl.Sb).put(0, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Te).put(0, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.I).put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Xe).put(0, new Integer[]{0,2,4,6,8});
		possibleStableValencies.get(ChemEl.Cs).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Ba).put(0, new Integer[]{2});
		possibleStableValencies.get(ChemEl.Tl).put(0, new Integer[]{1,3});
		possibleStableValencies.get(ChemEl.Pb).put(0, new Integer[]{2,4});
		possibleStableValencies.get(ChemEl.Bi).put(0, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Po).put(0, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.At).put(0, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Rn).put(0, new Integer[]{0,2,4,6,8});
		possibleStableValencies.get(ChemEl.Fr).put(0, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Ra).put(0, new Integer[]{2});

		possibleStableValencies.get(ChemEl.H).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Li).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Be).put(1, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Be).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.B).put(2, new Integer[]{1});
		possibleStableValencies.get(ChemEl.B).put(1, new Integer[]{2});
		possibleStableValencies.get(ChemEl.B).put(-1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.B).put(-2, new Integer[]{3});
		possibleStableValencies.get(ChemEl.C).put(2, new Integer[]{2});
		possibleStableValencies.get(ChemEl.C).put(1, new Integer[]{3});
		possibleStableValencies.get(ChemEl.C).put(-1, new Integer[]{3});
		possibleStableValencies.get(ChemEl.C).put(-2, new Integer[]{2});
		possibleStableValencies.get(ChemEl.N).put(2, new Integer[]{3});
		possibleStableValencies.get(ChemEl.N).put(1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.N).put(-1, new Integer[]{2});
		possibleStableValencies.get(ChemEl.N).put(-2, new Integer[]{1});
		possibleStableValencies.get(ChemEl.O).put(1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.O).put(1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.O).put(-1, new Integer[]{1});
		possibleStableValencies.get(ChemEl.O).put(-2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.F).put(2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.F).put(1, new Integer[]{2});
		possibleStableValencies.get(ChemEl.F).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Na).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Na).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Mg).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Al).put(3, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Al).put(2, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Al).put(1, new Integer[]{2});
		possibleStableValencies.get(ChemEl.Al).put(-1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.Al).put(-2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Si).put(2, new Integer[]{2});
		possibleStableValencies.get(ChemEl.Si).put(1, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Si).put(-1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Si).put(-2, new Integer[]{2});
		possibleStableValencies.get(ChemEl.P).put(2, new Integer[]{3});
		possibleStableValencies.get(ChemEl.P).put(1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.P).put(-1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.P).put(-2, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.S).put(2, new Integer[]{4});
		possibleStableValencies.get(ChemEl.S).put(1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.S).put(-1, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.S).put(-2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Cl).put(2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Cl).put(1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Cl).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.K).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.K).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ca).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ca).put(1, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Ga).put(3, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ga).put(2, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Ga).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ga).put(-1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.Ga).put(-2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Ge).put(4, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ge).put(1, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Ge).put(-1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Ge).put(-2, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.As).put(2, new Integer[]{3});
		possibleStableValencies.get(ChemEl.As).put(1, new Integer[]{4});
		possibleStableValencies.get(ChemEl.As).put(-1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.As).put(-2, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.As).put(-3, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Se).put(2, new Integer[]{4});
		possibleStableValencies.get(ChemEl.Se).put(1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Se).put(-1, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Se).put(-2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Br).put(2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Br).put(1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Br).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Rb).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Rb).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Sr).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Sr).put(1, new Integer[]{1});
		possibleStableValencies.get(ChemEl.In).put(3, new Integer[]{0});
		possibleStableValencies.get(ChemEl.In).put(2, new Integer[]{1});
		possibleStableValencies.get(ChemEl.In).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.In).put(-1, new Integer[]{2,4});
		possibleStableValencies.get(ChemEl.In).put(-2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Sn).put(4, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Sn).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Sn).put(1, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Sn).put(-1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Sn).put(-2, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Sb).put(3, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Sb).put(2, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Sb).put(1, new Integer[]{2,4});
		possibleStableValencies.get(ChemEl.Sb).put(-1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Sb).put(-2, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Te).put(2, new Integer[]{2,4});
		possibleStableValencies.get(ChemEl.Te).put(1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Te).put(-1, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.Te).put(-2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.I).put(2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.I).put(1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.I).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Cs).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Cs).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ba).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ba).put(1, new Integer[]{1});
		possibleStableValencies.get(ChemEl.Pb).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Pb).put(1, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Pb).put(-1, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.Pb).put(-2, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Bi).put(3, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Bi).put(2, new Integer[]{3});
		possibleStableValencies.get(ChemEl.Bi).put(1, new Integer[]{2,4});
		possibleStableValencies.get(ChemEl.Bi).put(-1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.Bi).put(-2, new Integer[]{1,3,5,7});
		possibleStableValencies.get(ChemEl.At).put(2, new Integer[]{3,5});
		possibleStableValencies.get(ChemEl.At).put(1, new Integer[]{2,4,6});
		possibleStableValencies.get(ChemEl.At).put(-1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Fr).put(1, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ra).put(2, new Integer[]{0});
		possibleStableValencies.get(ChemEl.Ra).put(1, new Integer[]{1});
	}

	/**
	 * Given a chemical element (e.g. Na) and charge (e.g. 1) returns the highest stable valency that OPSIN knows is possible
	 * If for the particular combination of chemical element and charge the highest stable valency is not known null is returned
	 * @param chemEl
	 * @param charge
	 * @return
	 */
	static Integer getMaximumValency(ChemEl chemEl, int charge) {
		Map<Integer, Integer[]> possibleStableValenciesForEl =  possibleStableValencies.get(chemEl);
		if (possibleStableValenciesForEl != null){
			Integer[] possibleStableValenciesForElAndCharge = possibleStableValenciesForEl.get(charge);
			if (possibleStableValenciesForElAndCharge != null){
				return possibleStableValenciesForElAndCharge[possibleStableValenciesForElAndCharge.length - 1];
			}
		}
		return null;
	}

	/**
	 * Checks whether the total incoming valency to an atom exceeds its expected valency
	 * outValency e.g. on radicals is taken into account
	 * @param a
	 * @return
	 */
	static boolean checkValency(Atom a) {
		int valency = a.getIncomingValency() + a.getOutValency();
		Integer maxVal;
		if (a.getLambdaConventionValency() != null){
			maxVal=a.getLambdaConventionValency() + a.getProtonsExplicitlyAddedOrRemoved();
		}
		else{
			ChemEl chemEl = a.getElement();
			int charge = a.getCharge();
			maxVal = getMaximumValency(chemEl, charge);
			if(maxVal == null) {
				return true;
			}
		}
		return valency <= maxVal;
	}

	/** Check whether valency is available on the atom to form a bond of the given order.
	 * spareValency and outValency are not taken into account.
	 * @param a atom you are interested in
	 * @param bondOrder order of bond required
     * @return
	 */
	static boolean checkValencyAvailableForBond(Atom a, int bondOrder) {
		int valency =a.getIncomingValency() + bondOrder;
		Integer maxVal;
		if (a.getLambdaConventionValency() != null){
			maxVal = a.getLambdaConventionValency() + a.getProtonsExplicitlyAddedOrRemoved();
		}
		else{
			ChemEl chemEl = a.getElement();
			int charge = a.getCharge();
			maxVal = getMaximumValency(chemEl, charge);
			if(maxVal == null) {
				return true;
			}
		}
		return valency <= maxVal;
	}

	/** Check whether changing to a heteroatom will result in valency being exceeded
	 * spareValency and outValency is taken into account
	 * @param a atom you are interested in
	 * @param heteroatom atom which will be replacing it
     * @return
	 */
	static boolean checkValencyAvailableForReplacementByHeteroatom(Atom a, Atom heteroatom) {
		int valency =a.getIncomingValency();
		valency +=a.hasSpareValency() ? 1 : 0;
		valency +=a.getOutValency();
		Integer maxValOfHeteroAtom = getMaximumValency(heteroatom.getElement(), heteroatom.getCharge());
        return maxValOfHeteroAtom == null || valency <= maxValOfHeteroAtom;
	}

	/**
	 * Returns the default valency of an element when uncharged or null if unknown
	 * @param chemlEl
	 * @return
	 */
	static Integer getDefaultValency(ChemEl chemlEl) {
		return expectedDefaultValency.get(chemlEl);
	}

	/**
	 * Returns the valency of an element in the HW system (useful for deciding whether something should have double bonds in a ring) or null if unknown
	 * Note that the HW system makes no claim about valency when the atom is charged
	 * @param chemEl
	 * @return
	 */
	static Integer getHWValency(ChemEl chemEl) {
		return valencyInHW.get(chemEl);
	}

	/**
	 * Returns the maximum valency of an element with a given charge or null if unknown
	 * @param chemEl
	 * @param charge
	 * @return
	 */
	static Integer[] getPossibleValencies(ChemEl chemEl, int charge) {
		Map<Integer, Integer[]> possibleStableValenciesForEl =  possibleStableValencies.get(chemEl);
		if (possibleStableValenciesForEl == null){
			return null;
		}
		return possibleStableValenciesForEl.get(charge);
	}
}
