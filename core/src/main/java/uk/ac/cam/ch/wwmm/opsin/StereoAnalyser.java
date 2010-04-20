package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

class StereoAnalyser {
	private static Map<String, Integer> elementToAtomicNumber = new HashMap<String, Integer>();

	/** Maps each atom to its currently assigned colour. Eventually all atoms in non identical environments will have different colours. Higher is higher priority*/
	private final Map<Atom, Integer> mappingToColour;

	/** Maps each atom to a list of of the colours of its neighbours*/
	private final Map<Atom, List<Integer>> atomNeighbourColours;

	/** The molecule upon which this StereoAnalyser is operating */
	private final Fragment molecule;

	private final SortByNeighbourColours sortByNeighbourColours;
	
	static{
		elementToAtomicNumber.put("H", 1);
		elementToAtomicNumber.put("He", 2);
		elementToAtomicNumber.put("Li", 3);
		elementToAtomicNumber.put("Be", 4);
		elementToAtomicNumber.put("B", 5);
		elementToAtomicNumber.put("C", 6);
		elementToAtomicNumber.put("N", 7);
		elementToAtomicNumber.put("O", 8);
		elementToAtomicNumber.put("F", 9);
		elementToAtomicNumber.put("Ne", 10);
		elementToAtomicNumber.put("Na", 11);
		elementToAtomicNumber.put("Mg", 12);
		elementToAtomicNumber.put("Al", 13);
		elementToAtomicNumber.put("Si", 14);
		elementToAtomicNumber.put("P", 15);
		elementToAtomicNumber.put("S", 16);
		elementToAtomicNumber.put("Cl", 17);
		elementToAtomicNumber.put("Ar", 18);
		elementToAtomicNumber.put("K", 19);
		elementToAtomicNumber.put("Ca", 20);
		elementToAtomicNumber.put("Sc", 21);
		elementToAtomicNumber.put("Ti", 22);
		elementToAtomicNumber.put("V", 23);
		elementToAtomicNumber.put("Cr", 24);
		elementToAtomicNumber.put("Mn", 25);
		elementToAtomicNumber.put("Fe", 26);
		elementToAtomicNumber.put("Co", 27);
		elementToAtomicNumber.put("Ni", 28);
		elementToAtomicNumber.put("Cu", 29);
		elementToAtomicNumber.put("Zn", 30);
		elementToAtomicNumber.put("Ga", 31);
		elementToAtomicNumber.put("Ge", 32);
		elementToAtomicNumber.put("As", 33);
		elementToAtomicNumber.put("Se", 34);
		elementToAtomicNumber.put("Br", 35);
		elementToAtomicNumber.put("Kr", 36);
		elementToAtomicNumber.put("Rb", 37);
		elementToAtomicNumber.put("Sr", 38);
		elementToAtomicNumber.put("Y", 39);
		elementToAtomicNumber.put("Zr", 40);
		elementToAtomicNumber.put("Nb", 41);
		elementToAtomicNumber.put("Mo", 42);
		elementToAtomicNumber.put("Tc", 43);
		elementToAtomicNumber.put("Ru", 44);
		elementToAtomicNumber.put("Rh", 45);
		elementToAtomicNumber.put("Pd", 46);
		elementToAtomicNumber.put("Ag", 47);
		elementToAtomicNumber.put("Cd", 48);
		elementToAtomicNumber.put("In", 49);
		elementToAtomicNumber.put("Sn", 50);
		elementToAtomicNumber.put("Sb", 51);
		elementToAtomicNumber.put("Te", 52);
		elementToAtomicNumber.put("I", 53);
		elementToAtomicNumber.put("Xe", 54);
		elementToAtomicNumber.put("Cs", 55);
		elementToAtomicNumber.put("Ba", 56);
		elementToAtomicNumber.put("La", 57);
		elementToAtomicNumber.put("Ce", 58);
		elementToAtomicNumber.put("Pr", 59);
		elementToAtomicNumber.put("Nd", 60);
		elementToAtomicNumber.put("Pm", 61);
		elementToAtomicNumber.put("Sm", 62);
		elementToAtomicNumber.put("Eu", 63);
		elementToAtomicNumber.put("Gd", 64);
		elementToAtomicNumber.put("Tb", 65);
		elementToAtomicNumber.put("Dy", 66);
		elementToAtomicNumber.put("Ho", 67);
		elementToAtomicNumber.put("Er", 68);
		elementToAtomicNumber.put("Tm", 69);
		elementToAtomicNumber.put("Yb", 70);
		elementToAtomicNumber.put("Lu", 71);
		elementToAtomicNumber.put("Hf", 72);
		elementToAtomicNumber.put("Ta", 73);
		elementToAtomicNumber.put("W", 74);
		elementToAtomicNumber.put("Re", 75);
		elementToAtomicNumber.put("Os", 76);
		elementToAtomicNumber.put("Ir", 77);
		elementToAtomicNumber.put("Pt", 78);
		elementToAtomicNumber.put("Au", 79);
		elementToAtomicNumber.put("Hg", 80);
		elementToAtomicNumber.put("Tl", 81);
		elementToAtomicNumber.put("Pb", 82);
		elementToAtomicNumber.put("Bi", 83);
		elementToAtomicNumber.put("Po", 84);
		elementToAtomicNumber.put("At", 85);
		elementToAtomicNumber.put("Rn", 86);
		elementToAtomicNumber.put("Fr", 87);
		elementToAtomicNumber.put("Ra", 88);
		elementToAtomicNumber.put("Ac", 89);
		elementToAtomicNumber.put("Th", 90);
		elementToAtomicNumber.put("Pa", 91);
		elementToAtomicNumber.put("U", 92);
		elementToAtomicNumber.put("Np", 93);
		elementToAtomicNumber.put("Pu", 94);
		elementToAtomicNumber.put("Am", 95);
		elementToAtomicNumber.put("Cm", 96);
		elementToAtomicNumber.put("Bk", 97);
		elementToAtomicNumber.put("Cf", 98);
		elementToAtomicNumber.put("Es", 99);
		elementToAtomicNumber.put("Fm", 100);
		elementToAtomicNumber.put("Md", 101);
		elementToAtomicNumber.put("No", 102);
		elementToAtomicNumber.put("Lr", 103);
		elementToAtomicNumber.put("Rf", 104);
		elementToAtomicNumber.put("Db", 105);
		elementToAtomicNumber.put("Sg", 106);
		elementToAtomicNumber.put("Bh", 107);
		elementToAtomicNumber.put("Hs", 108);
		elementToAtomicNumber.put("Mt", 109);
		elementToAtomicNumber.put("Ds", 110);
	}
	
	/**
	 * Holds information about a tetrahedral stereocentre
	 * @author dl387
	 *
	 */
	class StereoCentre{
		private final Atom stereoAtom;
		private final List<Atom> cipOrderedAtoms;
		/**
		 * Creates a stereocentre object from a tetrahedral stereocentre atom and a list of neighbouring atoms in CIP priority (lowest to highest)
		 * @param stereoAtom
		 * @param cipOrderedNeighbours
		 */
		StereoCentre(Atom stereoAtom, List<Atom> cipOrderedNeighbours) {
			this.stereoAtom = stereoAtom;
			this.cipOrderedAtoms = cipOrderedNeighbours;
		}
		public Atom getStereoAtom() {
			return stereoAtom;
		}
		public List<Atom> getCipOrderedAtoms() {
			return cipOrderedAtoms;
		}
	}
	
	/***
	 * Holds information about a double bond that can possess E/Z stereochemistry
	 * @author dl387
	 *
	 */
	class StereoBond{
		private final Bond bond;
		private final List<Atom> stereoAtoms;
		StereoBond(Bond bond, List<Atom> stereoAtoms) {
			this.bond = bond;
			this.stereoAtoms = stereoAtoms;
		}
		Bond getBond() {
			return bond;
		}
		
		/**
		 * Returns the following atoms:
		 * Highest CIP atom on one side
		 * atom in bond
		 * other atom in bond
		 * Highest CIP atom on other side
		 * @return
		 */
		List<Atom> getStereoAtoms() {
			return stereoAtoms;
		}
	}
	
	/**
	 * Sorts atoms by their colour, low to high
	 * Higher colour means higher CIP priority
	 * @author dl387
	 *
	 */
	private class SortByColour implements Comparator<Atom> {

	    public int compare(Atom a, Atom b){
	    	int colour1 = mappingToColour.get(a);
	    	int colour2 = mappingToColour.get(b);
	    	if (colour1 > colour2){
	    		return 1;
	    	}
	    	else if (colour1 < colour2){
	    		return -1;
	    	}
			return 0;
	    }
	}
	
	/**
	 * Sorts atoms by their atomic number, low to high
	 * @author dl387
	 *
	 */
	private static class SortAtomsByAtomicNumber implements Comparator<Atom> {

	    public int compare(Atom a, Atom b){
	    	int atomicNumber1 = elementToAtomicNumber.get(a.getElement());
	    	int atomicNumber2 = elementToAtomicNumber.get(b.getElement());
	    	if (atomicNumber1 > atomicNumber2){
	    		return 1;
	    	}
	    	else if (atomicNumber1 < atomicNumber2){
	    		return -1;
	    	}
			return 0;
	    }
	}
	
	/**
	 * Initially sorts on the atoms' colour and if these are the same then
	 * sorts based on the list of colours for neighbouring atoms 
	 * e.g. [1,2] > [1,1]  [1,1,3] > [2,2,2]  [1,1,3] > [3]  
	 * @author dl387
	 *
	 */
	private class SortByNeighbourColours implements Comparator<Atom> {
	    public int compare(Atom a, Atom b){
	    	int colour1 = mappingToColour.get(a);
	    	int colour2 = mappingToColour.get(b);
	    	if (colour1 > colour2){
	    		return 1;
	    	}
	    	else if (colour1 < colour2){
	    		return -1;
	    	}
	    	List<Integer> colours1 = atomNeighbourColours.get(a);
	    	List<Integer> colours2 = atomNeighbourColours.get(b);
	    	
	    	int colours1Size = colours1.size();
	    	int colours2Size = colours2.size();
	    	int differenceInSize = colours1Size - colours2Size;
	    	int maxCommonColourSize = colours1Size > colours2Size ? colours2Size : colours1Size;
	    	for (int i = 1; i <= maxCommonColourSize; i++) {
				int difference = colours1.get(colours1Size -i) - colours2.get(colours2Size -i);
				if (difference >0){
					return 1;
				}
				if (difference < 0){
					return -1;
				}
			}
	    	if (differenceInSize >0){
	    		return 1;
	    	}
	    	if (differenceInSize <0){
	    		return -1;
	    	}
			return 0;
	    }
	}
	/**
	 * Employs a derivative of the InChI algorithm to label which atoms are equivalent.
	 * These labels can then be used by the findStereo(Atoms/Bonds) functions to find features that
	 * can possess stereoChemistry
	 * @param molecule
	 * @throws StructureBuildingException
	 */
	StereoAnalyser(Fragment molecule) throws StructureBuildingException {
		this.molecule = molecule;
		sortByNeighbourColours = new SortByNeighbourColours();
		addGhostAtoms();
		List<Atom> atomList = molecule.getAtomList();
		mappingToColour = new HashMap<Atom, Integer>(atomList.size());
		atomNeighbourColours = new HashMap<Atom,List<Integer>>(atomList.size());
		Collections.sort(atomList, new SortAtomsByAtomicNumber());
		populateColoursByAtomicNumber(atomList);
		
		boolean changeFound = true;
		while(changeFound){
			for (Atom atom : atomList) {
				List<Integer> neighbourColours = findColourOfNeighbours(atom);
				atomNeighbourColours.put(atom, neighbourColours);
			}
			Collections.sort(atomList, sortByNeighbourColours);
			changeFound = populateColoursAndReportIfColoursWereChanged(atomList);
		}
		removeGhostAtoms();
	}

	/**
	 * Adds "ghost" atoms in accordance with the CIP rules for handling double bonds
	 * e.g. C=C --> C(G)=C(G) where ghost is a carbon with no hydrogen bonded to it
	 * @throws StructureBuildingException
	 */
	private void addGhostAtoms() throws StructureBuildingException {
		Set<Bond> bonds = molecule.getBondSet();
		int ghostIdCounter = -1;
		for (Bond bond : bonds) {
			int bondOrder = bond.getOrder();
			for (int i = bondOrder; i >1; i--) {
				Atom fromAtom =bond.getFromAtom();
				Atom toAtom =bond.getToAtom();

				Atom ghost1 = new Atom(ghostIdCounter--, fromAtom.getElement(), molecule);
				Bond b1 = new Bond(ghost1, toAtom, 1);
				toAtom.addBond(b1);
				ghost1.addBond(b1);
				molecule.addAtom(ghost1);
				
				Atom ghost2 = new Atom(ghostIdCounter--, toAtom.getElement(), molecule);
				Bond b2 = new Bond(ghost2, fromAtom, 1);
				fromAtom.addBond(b2);
				ghost2.addBond(b2);
				molecule.addAtom(ghost2);
			}
		}
	}

	/**
	 * Removes the ghost atoms added by addGhostAtoms
	 * @throws StructureBuildingException
	 */
	private void removeGhostAtoms() throws StructureBuildingException {
		List<Atom> atomList = molecule.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getID() < 0){
				Atom adjacentAtom = atom.getAtomNeighbours().get(0);
				adjacentAtom.removeBond(atom.getFirstBond());
				molecule.removeAtom(atom);
			}
		}
	}


	/**
	 * Takes a list of atoms sorted by atomicNumber
	 * and populates the mappingToColour map
	 * @param atomList
	 */
	private void populateColoursByAtomicNumber(List<Atom> atomList) {
		String lastAtomElement = atomList.get(0).getElement();
		List<Atom> atomsOfThisColour = new ArrayList<Atom>();
		int atomsSeen =0;
		for (Atom atom : atomList) {
			if (!atom.getElement().equals(lastAtomElement)){
				for (Atom a2 : atomsOfThisColour) {
					mappingToColour.put(a2, atomsSeen);
				}
				lastAtomElement = atom.getElement();
				atomsOfThisColour = new ArrayList<Atom>();
			}
			atomsOfThisColour.add(atom);
			atomsSeen++;
		}
		if (!atomsOfThisColour.isEmpty()){
			for (Atom a2 : atomsOfThisColour) {
				mappingToColour.put(a2, atomsSeen);
			}
		}
	}
	
	/**
	 * Takes a list of atoms sorted by colour/the colour of their neighbours
	 * and populates the mappingToColour map
	 * Returns whether mappingToColour was changed
	 * @param atomList
	 * @return boolean Whether mappingToColour was changed
	 */
	private boolean populateColoursAndReportIfColoursWereChanged(List<Atom> atomList) {
		Atom previousAtom = atomList.get(0);
		List<Atom> atomsOfThisColour = new ArrayList<Atom>();
		int atomsSeen =0;
		boolean changeFound = false;
		for (Atom atom : atomList) {
			if (sortByNeighbourColours.compare(previousAtom, atom)!=0){
				for (Atom atomOfThisColour : atomsOfThisColour) {
					if (!changeFound && atomsSeen != mappingToColour.get(atomOfThisColour)){
						changeFound =true;
					}
					mappingToColour.put(atomOfThisColour, atomsSeen);
				}
				previousAtom = atom;
				atomsOfThisColour = new ArrayList<Atom>();
			}
			atomsOfThisColour.add(atom);
			atomsSeen++;
		}
		if (!atomsOfThisColour.isEmpty()){
			for (Atom atomOfThisColour : atomsOfThisColour) {
				if (!changeFound && atomsSeen != mappingToColour.get(atomOfThisColour)){
					changeFound =true;
				}
				mappingToColour.put(atomOfThisColour, atomsSeen);
			}
		}
		return changeFound;
	}

	/**
	 * Produces a sorted (low to high) list of the colour of the atoms surrounding a given atom
	 * @param atom
	 * @return List<Integer> colourOfAdjacentAtoms
	 * @throws StructureBuildingException
	 */
	private List<Integer> findColourOfNeighbours(Atom atom) throws StructureBuildingException {	
		List<Integer> colourOfAdjacentAtoms = new ArrayList<Integer>();
		Set<Bond> bonds = atom.getBonds();
		for (Bond bond : bonds) {
			Atom otherAtom = bond.getFromAtom() == atom ? bond.getToAtom() : bond.getFromAtom();
			colourOfAdjacentAtoms.add(mappingToColour.get(otherAtom));
		} 
		Collections.sort(colourOfAdjacentAtoms);//sort such that this goes from low to high
		return colourOfAdjacentAtoms;
	}

	/**
	 * Retrieves a list of any tetrahedral stereoCentres
	 * Internally this is done by checking whether the "colour" of all neighbouring atoms of the tetrahedral atom are different
	 * @return List<StereoCentre>
	 * @throws StructureBuildingException
	 */
	List<StereoCentre> findStereoCentres() throws StructureBuildingException {
		List<Atom> atomList = molecule.getAtomList();
		List<StereoCentre> stereoCentres = new ArrayList<StereoCentre>();
		for (Atom atom : atomList) {
			List<Atom> neighbours = atom.getAtomNeighbours();
			if (isTetrahedral(atom)){
				int[] colours = new int[4];
				for (int i = neighbours.size() -1 ; i >=0; i--) {
					colours[i] = mappingToColour.get(neighbours.get(i));
				}
				
				boolean foundIdenticalNeighbour =false;
				for (int i = 0; i < 4; i++) {
					int cl = colours[i];
					for (int j = i +1; j < 4; j++) {
						if (cl == colours[j]){
							foundIdenticalNeighbour =true;
							break;
						}
					}
				}
				if (!foundIdenticalNeighbour){
					Collections.sort(neighbours, new SortByColour());
					if (neighbours.size()==3){//lone pair is the 4th. This is represented by the atom itself and is always the lowest priority
						neighbours.add(0, atom);
					}
					stereoCentres.add(new StereoCentre(atom, neighbours));
				}
			}
		}
		return stereoCentres;
	}
	
	/**
	 * Crudely determines whether an atom possesses tetrahedral geometry
	 * @param atom
	 * @return
	 * @throws StructureBuildingException
	 */
	private boolean isTetrahedral(Atom atom) throws StructureBuildingException {
		int neighbourCount = atom.getAtomNeighbours().size();
		String element = atom.getElement();
		if (neighbourCount == 4){
			if (element.equals("C") || element.equals("N")|| element.equals("P") || element.equals("S")
					|| element.equals("B")|| element.equals("Si") || element.equals("As") || element.equals("Se")){
				return true;
			}
		}
		else if (neighbourCount ==3){
			if ((element.equals("S") || element.equals("Se")) && atom.getIncomingValency()==4){//tetrahedral sulfur - 3 bonds and the lone pair
				return true;
			}
		}
		return false;
	}


	/**
	 *  Retrieves a list of any double bonds possessing the potential to have E/Z stereoChemistry
	 *  This is done internally by checking the two atoms attached to the ends of the double bond are different
	 *  As an exception nitrogen's lone pair is treated like a low priority group and so is allowed to only have 1 atom connected to it
	 * @return
	 * @throws StructureBuildingException
	 */
	List<StereoBond> findStereoBonds() throws StructureBuildingException {
		Set<Bond> bondSet =molecule.getBondSet();
		List<StereoBond> stereoBonds = new ArrayList<StereoBond>();
		for (Bond bond : bondSet) {
			if (bond.getOrder()==2){
				Atom a1 = bond.getFromAtom();
				List<Atom> neighbours1 = a1.getAtomNeighbours();
				neighbours1.remove(bond.getToAtom());
				if (neighbours1.size()==2 || (neighbours1.size()==1 && a1.getElement().equals("N") && a1.getIncomingValency()==3 && a1.getCharge()==0)){
					if (neighbours1.size()==2 && mappingToColour.get(neighbours1.get(0))== mappingToColour.get(neighbours1.get(1))){
						continue;
					}
					Atom a2 = bond.getToAtom();
					List<Atom> neighbours2 = a2.getAtomNeighbours();
					neighbours2.remove(bond.getFromAtom());
					if (neighbours2.size()==2 || (neighbours2.size()==1 && a2.getElement().equals("N") && a2.getIncomingValency()==3 && a2.getCharge()==0)){
						if (neighbours2.size()==2 && mappingToColour.get(neighbours2.get(0))== mappingToColour.get(neighbours2.get(1))){
							continue;
						}
						Collections.sort(neighbours1, new SortByColour());
						Collections.sort(neighbours2, new SortByColour());
						List<Atom> stereoAtoms = new ArrayList<Atom>();
						stereoAtoms.add(neighbours1.get(neighbours1.size()-1));//highest CIP adjacent to a1
						stereoAtoms.add(a1);
						stereoAtoms.add(a2);
						stereoAtoms.add(neighbours2.get(neighbours2.size()-1));//highest CIP adjacent to a2
						stereoBonds.add(new StereoBond(bond, stereoAtoms));
					}
				}
			}
		}
		return stereoBonds;
	}

}
