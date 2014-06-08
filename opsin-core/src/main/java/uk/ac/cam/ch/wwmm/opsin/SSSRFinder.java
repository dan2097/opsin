package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Class for finding SSSR
 * The algorithm employed does not work in some corner cases
 * 
 * @author pm286
 * @author dl387
 *
 */
class SSSRFinder {

	/** get set of smallest rings.
	 * In corner cases the list of rings returned will not be the SSSR
	 * @param frag 
	 * @return list of rings
	 */
	static List<Ring> getSetOfSmallestRings(Fragment frag){
		List<Atom> atomSet = frag.getAtomList();
				
		List<Ring> ringList = getRings(atomSet);
		
		if (ringList.size() > 1) {
			boolean change = true;
			while (change) {
				for (int i = 0; i < ringList.size(); i++) {
					Ring ring = ringList.get(i);
					change = reduceRingSizes(ring, ringList);
				}
			}
		}		
		return ringList;
	}

	/** get list of rings.
	 * not necessarily SSSR
	 * @param atomSet
	 * @return list of rings
	 */
	private static List<Ring> getRings(List<Atom> atomSet ){
		List<Ring> ringList = new ArrayList<Ring>();
		Set<Atom> usedAtoms = new HashSet<Atom>();
		
		Atom root = atomSet.get(0); 
		Atom parentAtom = null;
		Map<Atom, Atom> atomToParentMap = new HashMap<Atom, Atom>();
		Set<Bond> linkBondSet = new LinkedHashSet<Bond>(); 
		
		expand(root, parentAtom, usedAtoms, atomToParentMap, linkBondSet);
		
		for (Bond bond : linkBondSet) {
			Ring ring = getRing(bond, atomToParentMap);
			ringList.add(ring);
		}
		
		return ringList;
	}
	
	private static Ring getRing(Bond bond, Map<Atom, Atom> atomToParentMap){ 
		Atom atomFrom =  bond.getFromAtom() ;
		Atom atomTo = bond.getToAtom(); 
		List<Bond>  bondSet0 = getAncestors1(atomFrom, atomToParentMap);
		List<Bond>  bondSet1 = getAncestors1(atomTo, atomToParentMap);
		List<Bond>  mergedBondSet = symmetricDifference (bondSet0, bondSet1); 
		
		mergedBondSet.add(bond);
		return new Ring(mergedBondSet);
	}
	
	private static List<Bond> getAncestors1(Atom atom, Map<Atom, Atom> atomToParentMap){
		List<Bond> newBondSet = new ArrayList<Bond>();
		while (true) {
			Atom atom1 = atomToParentMap.get(atom);
			if (atom1 == null) {
				break;
			}
			Bond bond = atom.getBondToAtom(atom1);
			if (newBondSet.contains(bond)) {
				break;
			}
			newBondSet.add(bond);
			atom = atom1;
		}
		return newBondSet;
	}
	
	private static void expand(Atom atom, Atom parentAtom, 
			Set<Atom> usedAtoms, Map<Atom, Atom> atomToParentMap,
			Set<Bond> linkBondSet){
		
		usedAtoms.add(atom);
		atomToParentMap.put(atom, parentAtom);
		List<Atom> ligandAtomList = atom.getAtomNeighbours();
		
		for (Atom ligandAtom : ligandAtomList) {
			if (ligandAtom.equals(parentAtom)) {
				// skip existing bond
			} else if (usedAtoms.contains(ligandAtom)) {
				Bond linkBond = atom.getBondToAtom(ligandAtom);
				linkBondSet.add(linkBond);
				// already treated
			} else {
				expand(ligandAtom, atom, usedAtoms, atomToParentMap, linkBondSet);
			}
		}
	}

	private static boolean reduceRingSizes(Ring ring, List<Ring> newList){
		boolean change = false;
		for (int i = 0; i < newList.size(); i++) {
			Ring target = newList.get(i);
			if (target == ring) {
				continue;
			}
			
			List<Bond> newBondSet = symmetricDifference ( target.getBondList(), ring.getBondList() ) ;
			if (newBondSet.size() < target.size()) {
				Ring newRing = new Ring(newBondSet);
				newList.set(i, newRing);
				change = true;
			}
		}
		return change;
	}

	
	private static List<Bond> symmetricDifference(List<Bond> bondSet1, List<Bond> bondSet2) {
		List<Bond> newBondSet = new ArrayList<Bond>();

		for (Bond bond1 : bondSet1) {
			if (!bondSet2.contains(bond1)) {
				newBondSet.add(bond1);
			}
		}
		for (Bond bond2 : bondSet2) {
			if (!bondSet1.contains(bond2)) {
				newBondSet.add(bond2);
			}
		}
		return newBondSet;
	}

}
