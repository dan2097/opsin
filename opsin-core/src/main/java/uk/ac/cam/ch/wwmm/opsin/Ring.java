package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

/**
 * Class representing a single ring (i.e. NOT a fused ring which is formed from multiple rings)
 * @author dl387
 *
 */
class Ring {
	private List<Atom> atomList = new ArrayList<Atom>();
	private List<Bond> bondList;
	private List<Atom> cyclicAtomList;
	private List<Bond> cyclicBondList;
	private List<Ring> neighbours = new ArrayList<Ring>();

	private int nFusedBonds = 0;
	
	Ring(List<Bond> bondList) throws StructureBuildingException{
		if (bondList==null || bondList.size()<=0){
			throw new StructureBuildingException("Bond list is empty");
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
			throw new StructureBuildingException("atomList and bondList different sizes. Ring(bond)");
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

		for (Bond bond : bondList) {
			if (bond.getFusedRings().size()>0) {
				bonds.add(bond);
			}
		}
		return bonds;
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
	void makeCyclicLists(Bond stBond, Atom stAtom){
		if (cyclicBondList==null){
			cyclicBondList = new ArrayList<Bond>();
			cyclicAtomList = new ArrayList<Atom>();

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
		String st = "";
		List<Atom> toPrint;

		if (cyclicAtomList != null) {
			toPrint = cyclicAtomList;
		}
		else {
			toPrint = atomList;
		}
		for (Atom atom : toPrint){
			st+=atom.getID() + " ";
		}
		return st;
	}

}
