package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

/**
 * Class representing a single ring (i.e. NOT a fused ring which is formed from multiple rings)
 * @author dl387
 *
 */
class Ring {
	private List<Atom>  atomSet;
	private List<Bond>  bondSet;
	private List<Atom>  cyclicAtomSet;
	private List<Bond>  cyclicBondSet;
	private List<Ring>  neighbours = new ArrayList<Ring>();

	private int nFusedBonds = 0;

	Ring(List<Atom> atomSet) throws StructureBuildingException {
		this.atomSet = atomSet;
		List<Bond> orderedBonds = new ArrayList<Bond>();
		CyclicAtomList ringIteratorFromAtom0 = new CyclicAtomList(atomSet);
		CyclicAtomList ringIteratorFromAtom1 = new CyclicAtomList(atomSet, 0);
		for (int i = 0; i < atomSet.size(); i++) {
			Atom a1 = ringIteratorFromAtom0.getNext();
			Atom a2 = ringIteratorFromAtom1.getNext();
			orderedBonds.add(a1.getFrag().findBondOrThrow(a1, a2));
		}
		this.bondSet = orderedBonds;
	}


	List<Bond>  getBondSet() {
		return bondSet;
	}

	List<Atom>  getAtomSet() {
		return atomSet;
	}

	/**
	 * Number of ring atoms/bonds
	 * @return
	 */
	int size() {
		return atomSet.size();
	}

	int getNumberOfFusedBonds() {
		return nFusedBonds;
	}

	void incrementNumberOfFusedBonds() {
		this.nFusedBonds++;
	}

	/**
	 * Return bonds utilised in multiple rings
	 * @return List<Bond>
	 */
	List<Bond> getFusedBonds(){
		ArrayList<Bond> bonds = new ArrayList<Bond>();

		for (Bond bond : bondSet) {
			if (bond.getFusedRings().size()>0) bonds.add(bond);
		}
		return bonds;
	}

	int getBondNumber(Bond bond){
		return cyclicBondSet.indexOf(bond);
	}

	List<Bond> getCyclicBondSet(){
		return cyclicBondSet;
	}

	List<Atom> getCyclicAtomSet(){
		return cyclicAtomSet;
	}

	List<Ring> getNeighbours() {
		return neighbours;
	}

	void addNeighbour(Ring ring) {
		this.neighbours.add(ring);
	}

	/**
	 * Stores atoms and bonds in the order defined by atom and bond
	 * @param stBond - the first bond
	 * @param stAtom - the atom defining in which direction to go
	 */
	void makeCyclicSets(Bond stBond, Atom stAtom){
		if (cyclicBondSet==null){
			cyclicBondSet = new ArrayList<Bond>();
			cyclicAtomSet = new ArrayList<Atom>();

			Atom atom = stAtom;
			cyclicBondSet.add(stBond);
			cyclicAtomSet.add(atom);

			for (int i=0; i<size()-1; i++){
				for(Bond bond2 : bondSet){
					if (cyclicBondSet.contains(bond2)) continue;
					if (bond2.getFromAtom() == atom){
						cyclicBondSet.add(bond2);
						atom = bond2.getToAtom();
						cyclicAtomSet.add(atom);
					}
					else if (bond2.getToAtom() == atom){
						cyclicBondSet.add(bond2);
						atom = bond2.getFromAtom();
						cyclicAtomSet.add(atom);
					}
				}
			}
		}
	}

	public String toString() {
		String st = "";
		List<Atom> toPrint;

		if (cyclicAtomSet != null) {
			toPrint = cyclicAtomSet;
		}
		else {
			toPrint = atomSet;
		}
		for (Atom atom : toPrint){
			st+=atom.getID() + " ";
		}
		return st;
	}

}
