package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

/**
 * Class representing a single ring (i.e. NOT a fused ring which is formed from multiple rings)
 * @author dl387
 *
 */
class Ring {
	private List<Atom> atomSet = new ArrayList<Atom>();
	private List<Bond> bondSet;
	private List<Atom> cyclicAtomSet;
	private List<Bond> cyclicBondSet;
	private List<Ring> neighbours = new ArrayList<Ring>();

	private int nFusedBonds = 0;
	
	Ring(List<Bond> bondSet) throws StructureBuildingException{
		if (bondSet==null || bondSet.size()<=0) throw new StructureBuildingException("Bond set is empty");
		this.bondSet = bondSet;

		for(Bond bond: bondSet){
			Atom atom1 = bond.getFromAtom();				
			if (!atomSet.contains(atom1)) {
				atomSet.add(atom1);
			}
				
			Atom atom2 = bond.getToAtom();
			if (!atomSet.contains(atom2)) {
				atomSet.add(atom2);
			}
		}
		
		if (atomSet.size() != bondSet.size()) {
			throw new StructureBuildingException("atomSet and bondset different sizes. Ring(bond)");
		}
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
		nFusedBonds++;
	}

	/**
	 * Return bonds utilised in multiple rings
	 * @return List<Bond>
	 */
	List<Bond> getFusedBonds(){
		List<Bond> bonds = new ArrayList<Bond>();

		for (Bond bond : bondSet) {
			if (bond.getFusedRings().size()>0) {
				bonds.add(bond);
			}
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
		neighbours.add(ring);
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
					if (cyclicBondSet.contains(bond2)){
						continue;
					}
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
