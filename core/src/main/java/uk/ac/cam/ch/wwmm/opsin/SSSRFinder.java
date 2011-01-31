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
     * @throws StructureBuildingException
	 */
	static List<Ring> getSetOfSmallestRings(Fragment frag) throws StructureBuildingException {
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
		else throw new StructureBuildingException("Ring perception system found less than 2 rings within input fragment!");
				
		return ringList;
	}

	/** get list of rings.
	 * not necessarily SSSR
	 * @param atomSet
     * @return list of rings
	 * @throws StructureBuildingException 
	 */
	private static List<Ring> getRings(List<Atom> atomSet ) throws StructureBuildingException {
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
	
	private static Ring getRing(Bond bond, Map<Atom, Atom> atomToParentMap) throws StructureBuildingException { 
		Atom atomFrom =  bond.getFromAtom() ;
		Atom atomTo = bond.getToAtom(); 
		List<Atom>  atomSet0 = getAncestors(atomFrom, atomToParentMap);
		List<Atom>  atomSet1 = getAncestors(atomTo, atomToParentMap);
		List<Bond>  bondSet0 = getAncestors1(atomFrom, atomToParentMap, atomSet1);
		List<Bond>  bondSet1 = getAncestors1(atomTo, atomToParentMap, atomSet0);
		List<Bond>  mergedBondSet = symmetricDifference (bondSet0, bondSet1); 
		
		mergedBondSet.add(bond);
		Ring ring = new Ring(mergedBondSet);
		return ring;
	}
	
	
	private static List<Atom>  getAncestors(Atom atom, Map<Atom, Atom> atomToParentMap) {
		List<Atom> newAtomSet = new ArrayList<Atom> ();

		atom = atomToParentMap.get(atom);
		if (atom != null && !newAtomSet.contains(atom)) {
			newAtomSet.add(atom);
		}

		return newAtomSet;
	}
	
	private static List<Bond>  getAncestors1(Atom atom, Map<Atom, Atom> atomToParentMap, List<Atom> atomSet) throws StructureBuildingException {
		Fragment molecule =atom.getFrag();
		List<Bond> newBondSet = new ArrayList<Bond>();
		while (true) {
			Atom atom1 = atomToParentMap.get(atom);
			if (atom1 == null) {
				break;
			}
			Bond bond = molecule.findBondOrThrow(atom, atom1);
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
			Set<Bond> linkBondSet) throws StructureBuildingException {
		
		usedAtoms.add(atom);
		atomToParentMap.put(atom, parentAtom);
		List<Atom> ligandAtomList = atom.getAtomNeighbours();
		Fragment fragment = atom.getFrag();
				
		for (int i = 0; i < ligandAtomList.size(); i++) {
			Atom ligandAtom = ligandAtomList.get(i);
			if (ligandAtom.equals(parentAtom)) {
				// skip existing bond
			} else if (usedAtoms.contains(ligandAtom)) {
				Bond linkBond = fragment.findBondOrThrow(atom, ligandAtom);
				linkBondSet.add(linkBond);
				// already treated
			} else {
				expand(ligandAtom, atom, usedAtoms, atomToParentMap, linkBondSet);
			}
		}
	}
	
	
	/**
	 * @param ring
	 * @param newList
     * @throws StructureBuildingException
     * @return
	 */
	private static boolean reduceRingSizes(Ring ring, List<Ring> newList) throws StructureBuildingException {
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
	        
        for (int i = 0; i < bondSet1.size(); i++) {
            if (!bondSet2.contains(bondSet1.get(i))) {
                newBondSet.add(bondSet1.get(i));
            }
        }
        for (int i = 0; i < bondSet2.size(); i++) {
            Bond bond = bondSet2.get(i);
            if (!bondSet1.contains(bond)) {
                newBondSet.add(bond);
            }
        }

	    return newBondSet;
	 }

}
