package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class representing a single ring (i.e. NOT a fused ring which is formed from multiple rings)
 * @author dl387
 *
 */
class Ring {
	private final List<Atom> atomList = new ArrayList<>();
	private final List<Bond> bondList;
	private final Map<Bond, Ring> bondToNeighbourRings = new LinkedHashMap<>();
	
	private List<Atom> cyclicAtomList;
	private List<Bond> cyclicBondList;
	
	Ring(List<Bond> bondList){
		if (bondList==null || bondList.size()==0){
			throw new IllegalArgumentException("Bond list is empty");
		}
		this.bondList = bondList;

		for(Bond bond: bondList){
			Atom atom1 = bond.getFromAtom();				
			if (!atomList.contains(atom1)) {
				atomList.add(atom1);
			}
				
			Atom atom2 = bond.getToAtom();
			if (!atomList.contains(atom2)) {
				atomList.add(atom2);
			}
		}
		
		if (atomList.size() != bondList.size()) {
			throw new RuntimeException("atomList and bondList different sizes. Ring(bond)");
		}
	}


	List<Bond>  getBondList() {
		return bondList;
	}

	List<Atom>  getAtomList() {
		return atomList;
	}

	/**
	 * Number of ring atoms/bonds
	 * @return
	 */
	int size() {
		return atomList.size();
	}

	int getNumberOfFusedBonds() {
		return bondToNeighbourRings.size();
	}

	/**
	 * Return bonds utilised in multiple rings
	 * @return List<Bond>
	 */
	List<Bond> getFusedBonds(){
		return new ArrayList<>(bondToNeighbourRings.keySet());
	}

	int getBondIndex(Bond bond){
		return cyclicBondList.indexOf(bond);
	}

	List<Bond> getCyclicBondList(){
		return cyclicBondList;
	}

	List<Atom> getCyclicAtomList(){
		return cyclicAtomList;
	}

	List<Ring> getNeighbours() {
		return new ArrayList<>(bondToNeighbourRings.values());
	}
	
	Ring getNeighbourOfFusedBond(Bond fusedBond) {
		return bondToNeighbourRings.get(fusedBond);
	}

	void addNeighbour(Bond bond, Ring ring) {
		if (this == ring) {
			throw new IllegalArgumentException("Ring can't be a neighbour of itself");
		}
		bondToNeighbourRings.put(bond, ring);
	}

	/**
	 * Stores atoms and bonds in the order defined by atom and bond
	 * @param stBond - the first bond
	 * @param stAtom - the atom defining in which direction to go
	 */
	void makeCyclicLists(Bond stBond, Atom stAtom){
		if (cyclicBondList==null){
			cyclicBondList = new ArrayList<>();
			cyclicAtomList = new ArrayList<>();

			Atom atom = stAtom;
			cyclicBondList.add(stBond);
			cyclicAtomList.add(atom);

			for (int i=0; i<size()-1; i++){
				for(Bond bond2 : bondList){
					if (cyclicBondList.contains(bond2)){
						continue;
					}
					if (bond2.getFromAtom() == atom){
						cyclicBondList.add(bond2);
						atom = bond2.getToAtom();
						cyclicAtomList.add(atom);
					}
					else if (bond2.getToAtom() == atom){
						cyclicBondList.add(bond2);
						atom = bond2.getFromAtom();
						cyclicAtomList.add(atom);
					}
				}
			}
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Atom atom : atomList){
			sb.append(atom.getID());
			sb.append(" ");
		}
		return sb.toString();
	}

}
