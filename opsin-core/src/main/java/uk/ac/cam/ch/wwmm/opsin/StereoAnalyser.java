package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

/**
 * Identifies stereocentres and determines the CIP order of connected atoms
 * @author dl387
 *
 */
class StereoAnalyser {
	/** Maps each atom to its currently assigned colour. Eventually all atoms in non identical environments will have different colours. Higher is higher priority*/
	private final Map<Atom, Integer> mappingToColour;

	/** Maps each atom to a list of of the colours of its neighbours*/
	private final Map<Atom, List<Integer>> atomNeighbourColours;

	/** The molecule upon which this StereoAnalyser is operating */
	private final Fragment molecule;

	private final AtomColourThenNeighbouringColoursComparator atomColourThenNeighbouringColoursComparator;
	private final AtomicNumberComparator atomicNumberComparator;
	private final AtomListCIPComparator atomListCIPComparator;
	private final ListOfAtomListsCIPComparator listOfAtomListsCIPComparator;
	
	/** The id of the next ghost to create*/
	private int ghostIdCounter = -1;
	
	//phantom atoms are not added as I believe that the results of the program will still be the same even in their absence as everything beats a phantom and comparing phantoms to phantoms achieves nothing
	
	/**
	 * Holds information about a tetrahedral stereocentre
	 * @author dl387
	 *
	 */
	class StereoCentre{
		private final Atom stereoAtom;
		/**
		 * Creates a stereocentre object from a tetrahedral stereocentre atom
		 * @param stereoAtom
		 */
		StereoCentre(Atom stereoAtom) {
			this.stereoAtom = stereoAtom;
		}
		Atom getStereoAtom() {
			return stereoAtom;
		}
		List<Atom> getCipOrderedAtoms() throws StructureBuildingException {
			List<Atom> cipOrderedAtoms = getNeighbouringAtomsInCIPOrder(stereoAtom);
			if (cipOrderedAtoms.size()==3){//lone pair is the 4th. This is represented by the atom itself and is always the lowest priority
				cipOrderedAtoms.add(0, stereoAtom);
			}
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
		StereoBond(Bond bond) {
			this.bond = bond;
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
		 * @throws StructureBuildingException 
		 */
		List<Atom> getOrderedStereoAtoms() throws StructureBuildingException {
			Atom a1 = bond.getFromAtom();
			Atom a2 = bond.getToAtom();
			List<Atom> cipOrderedNeighbours1 = getNeighbouringAtomsInCIPOrder(a1);
			List<Atom> cipOrderedNeighbours2 = getNeighbouringAtomsInCIPOrder(a2);
			cipOrderedNeighbours1.remove(a2);
			cipOrderedNeighbours2.remove(a1);
			List<Atom> stereoAtoms = new ArrayList<Atom>();
			stereoAtoms.add(cipOrderedNeighbours1.get(cipOrderedNeighbours1.size()-1));//highest CIP adjacent to a1
			stereoAtoms.add(a1);
			stereoAtoms.add(a2);
			stereoAtoms.add(cipOrderedNeighbours2.get(cipOrderedNeighbours2.size()-1));//highest CIP adjacent to a2
			return stereoAtoms;
		}
	}
	
	/**
	 * Sorts atoms by their atomic number, low to high
	 * @author dl387
	 *
	 */
	private static class AtomicNumberComparator implements Comparator<Atom> {

	    public int compare(Atom a, Atom b){
	    	int atomicNumber1 = AtomProperties.elementToAtomicNumber.get(a.getElement());
	    	int atomicNumber2 = AtomProperties.elementToAtomicNumber.get(b.getElement());
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
	private class AtomColourThenNeighbouringColoursComparator implements Comparator<Atom> {
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
		atomColourThenNeighbouringColoursComparator = new AtomColourThenNeighbouringColoursComparator();
		atomicNumberComparator = new AtomicNumberComparator();
		atomListCIPComparator = new AtomListCIPComparator();
		listOfAtomListsCIPComparator = new ListOfAtomListsCIPComparator();
		addGhostAtoms();
		List<Atom> atomList = molecule.getAtomList();
		mappingToColour = new HashMap<Atom, Integer>(atomList.size());
		atomNeighbourColours = new HashMap<Atom,List<Integer>>(atomList.size());
		Collections.sort(atomList, atomicNumberComparator);
		populateColoursByAtomicNumber(atomList);
		
		boolean changeFound = true;
		while(changeFound){
			for (Atom atom : atomList) {
				List<Integer> neighbourColours = findColourOfNeighbours(atom);
				atomNeighbourColours.put(atom, neighbourColours);
			}
			Collections.sort(atomList, atomColourThenNeighbouringColoursComparator);
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
	 */
	private void removeGhostAtoms() {
		List<Atom> atomList = molecule.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getID() < 0){
				Atom adjacentAtom = atom.getAtomNeighbours().get(0);
				adjacentAtom.removeBond(atom.getFirstBond());
				molecule.removeAtom(atom);
			}
		}
		ghostIdCounter = -1;
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
			if (atomColourThenNeighbouringColoursComparator.compare(previousAtom, atom)!=0){
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
	 */
	private List<Integer> findColourOfNeighbours(Atom atom) {	
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
	 */
	List<StereoCentre> findStereoCentres(){
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
					stereoCentres.add(new StereoCentre(atom));
				}
			}
		}
		return stereoCentres;
	}
	
	/**
	 * Crudely determines whether an atom possesses tetrahedral geometry
	 * @param atom
	 * @return
	 */
	static boolean isTetrahedral(Atom atom) {
		int neighbourCount = atom.getAtomNeighbours().size();
		String element = atom.getElement();
		if (neighbourCount == 4){
			if (element.equals("B") || element.equals("C") || element.equals("Si") || element.equals("Ge") ||
					element.equals("N")|| element.equals("P") || element.equals("S") ||
					  element.equals("As") || element.equals("Se")){
				return true;
			}
		}
		else if (neighbourCount ==3){
			if ((element.equals("S") || element.equals("Se")) && (atom.getIncomingValency()==4 || (atom.getCharge() ==1 && atom.getIncomingValency()==3))){
				//tetrahedral sulfur/selenium - 3 bonds and the lone pair
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
	 */
	List<StereoBond> findStereoBonds() {
		Set<Bond> bondSet =molecule.getBondSet();
		List<StereoBond> stereoBonds = new ArrayList<StereoBond>();
		for (Bond bond : bondSet) {
			if (bond.getOrder()==2){
				Atom a1 = bond.getFromAtom();
				List<Atom> neighbours1 =  a1.getAtomNeighbours();
				neighbours1.remove(bond.getToAtom());
				if (neighbours1.size()==2 || (neighbours1.size()==1 && a1.getElement().equals("N") && a1.getIncomingValency()==3 && a1.getCharge()==0)){
					if (neighbours1.size()==2 && mappingToColour.get(neighbours1.get(0)).equals(mappingToColour.get(neighbours1.get(1)))){
						continue;
					}
					Atom a2 = bond.getToAtom();
					List<Atom> neighbours2 = a2.getAtomNeighbours();
					neighbours2.remove(bond.getFromAtom());
					if (neighbours2.size()==2 || (neighbours2.size()==1 && a2.getElement().equals("N") && a2.getIncomingValency()==3 && a2.getCharge()==0)){
						if (neighbours2.size()==2 && mappingToColour.get(neighbours2.get(0)).equals(mappingToColour.get(neighbours2.get(1)))){
							continue;
						}
						stereoBonds.add(new StereoBond(bond));
					}
				}
			}
		}
		return stereoBonds;
	}
	
	/**
	 * Returns the given atoms neighbours in CIP order from lowest priority to highest priority
	 * @param chiralAtom
	 * @return
	 * @throws StructureBuildingException
	 */
	List<Atom> getNeighbouringAtomsInCIPOrder(Atom chiralAtom) throws StructureBuildingException{
		List<Atom> neighbours = chiralAtom.getAtomNeighbours();
		addGhostAtomsForCIPAssignment(chiralAtom);
		Collections.sort(neighbours, new SortByCIPOrder(chiralAtom));
		removeGhostAtoms();
		return neighbours;
	}
	
	/**
	 * Adds "ghost" atoms in accordance with the CIP rules for handling double bonds
	 * e.g. C=C --> C(G)=C(G) where ghost is a carbon with no hydrogen bonded to it
	 * Higher order bonds connected to the chiral atom are not converted in accordance with P-91.1.4.2.4 (IUPAC 2004 guidelines)
	 * @param chiralAtom Higher order bonds connected to this atom are not touched
	 * @throws StructureBuildingException
	 */
	private void addGhostAtomsForCIPAssignment(Atom chiralAtom) throws StructureBuildingException {
		Set<Bond> bonds = molecule.getBondSet();
		for (Bond bond : bonds) {
			int bondOrder = bond.getOrder();
			for (int i = bondOrder; i >1; i--) {
				Atom fromAtom =bond.getFromAtom();
				Atom toAtom =bond.getToAtom();
				if (!fromAtom.equals(chiralAtom) && !toAtom.equals(chiralAtom)){//P-91.1.4.2.4
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
	}
	
	/**
	 * Sorts atoms by their CIP order, low to high
	 * @author dl387
	 *
	 */
	private class SortByCIPOrder implements Comparator<Atom> {
		private final Atom chiralAtom;
	    SortByCIPOrder(Atom chiralAtom) {
			this.chiralAtom = chiralAtom;
		}

		public int compare(Atom a, Atom b){
			int atomicNumber1 = AtomProperties.elementToAtomicNumber.get(a.getElement());
			int atomicNumber2 = AtomProperties.elementToAtomicNumber.get(b.getElement());
	    	if (atomicNumber1 > atomicNumber2){
	    		return 1;
	    	}
	    	else if (atomicNumber1 < atomicNumber2){
	    		return -1;
	    	}
	    	int currentGhostId = ghostIdCounter; 
			CipState startingState = prepareInitialCIPState(a, b);
	    	Queue<CipState> cipStateQueue = new LinkedList<CipState>();
	    	cipStateQueue.add(startingState);
	    	/* Go through CIP states in a breadth-first manner:
	    	 * Neighbours of the given atom/s (if multiple atoms this is because so far the two paths leading to them have been equivalent) are evaluated for both a and b
	    	 * Neighbours are sorted by CIP priority
	    	 * Comparisons performed between neighbours of a and neighbours of b (will break if compare!=0)
	    	 * Degenerate neighbours grouped together
	    	 * CIP state formed for each list of neighbours and added to queue in order of priority
	    	 *
	    	 */
	    	int compare =0;
	    	while(!cipStateQueue.isEmpty()){
	    		CipState currentState = cipStateQueue.remove();
	    		compare = compareAtNextLevel(currentState, cipStateQueue);
	    		if (compare!=0){
	    			break;
	    		}
	    	}
	    	removeAnyGhostAtomsAddedAndCorrectGhostIdCounter(currentGhostId);
	    	return compare;
	    }

		private CipState prepareInitialCIPState(Atom a, Atom b) {	
			List<Atom> nextAtoms1 = new ArrayList<Atom>();
			nextAtoms1.add(a);
			List<Atom> nextAtoms2 = new ArrayList<Atom>();
			nextAtoms2.add(b);
			
			List<Atom> previousAtoms1 = new ArrayList<Atom>();
			previousAtoms1.add(chiralAtom);
			List<Atom> previousAtoms2 = new ArrayList<Atom>();
			previousAtoms2.add(chiralAtom);
			
			Map<Atom,Set<Atom>> previousAtomToVisitedAtoms1 = new HashMap<Atom, Set<Atom>>();
			Set<Atom> atomsVisted = new HashSet<Atom>();
			atomsVisted.add(chiralAtom);
			previousAtomToVisitedAtoms1.put(chiralAtom, atomsVisted);
			Map<Atom,Set<Atom>> previousAtomToVisitedAtoms2 = new HashMap<Atom, Set<Atom>>();
			previousAtomToVisitedAtoms2.put(chiralAtom, new HashSet<Atom>(atomsVisted));
	    	CipState startingState = new CipState(previousAtomToVisitedAtoms1, previousAtomToVisitedAtoms2, previousAtoms1, previousAtoms2, nextAtoms1, nextAtoms2);
			return startingState;
		}
	}
	
	/**
	 * Holds information about what atoms to try next next, how those atoms are reached (to prevent immediate back tracking and to detect cycles)
	 * @author dl387
	 *
	 */
	private static class CipState{
		CipState(Map<Atom, Set<Atom>> previousAtomToVisitedAtoms1, Map<Atom, Set<Atom>> previousAtomToVisitedAtoms2, 
				List<Atom> previousAtoms1, List<Atom> previousAtoms2, List<Atom> nextAtoms1, List<Atom> nextAtoms2) {
			this.previousAtomToVisitedAtoms1 = previousAtomToVisitedAtoms1;
			this.previousAtomToVisitedAtoms2 = previousAtomToVisitedAtoms2;
			this.previousAtoms1 = previousAtoms1;
			this.previousAtoms2 = previousAtoms2;
			this.nextAtoms1 = nextAtoms1;
			this.nextAtoms2 = nextAtoms2;
		}
		final Map<Atom,Set<Atom>> previousAtomToVisitedAtoms1;
		final Map<Atom,Set<Atom>> previousAtomToVisitedAtoms2;
		final List<Atom> previousAtoms1;
		final List<Atom> previousAtoms2;
		final List<Atom> nextAtoms1;
		final List<Atom> nextAtoms2;
	}
	
	/**
	 * Compares the neighbours of the atoms specified in nextAtom1/2 in cipstate.
	 * Returns the result of the comparison between these neighbours
	 * If the comparison returned 0 adds new cipstates to the queue
	 * @param cipState
	 * @param queue
	 * @return
	 */
	private int compareAtNextLevel(CipState cipState, Queue<CipState> queue) {
		Map<Atom, Atom> currentToPrevious1 = new HashMap<Atom, Atom>();
		List<List<List<Atom>>> newNeighbours1 = getNextLevelNeighbours(cipState.previousAtomToVisitedAtoms1, cipState.previousAtoms1, cipState.nextAtoms1, currentToPrevious1);
		Map<Atom, Atom> currentToPrevious2 = new HashMap<Atom, Atom>();
		List<List<List<Atom>>> newNeighbours2 = getNextLevelNeighbours(cipState.previousAtomToVisitedAtoms2, cipState.previousAtoms2, cipState.nextAtoms2, currentToPrevious2);

		int compare = compareNeighboursByCIPpriorityRules(newNeighbours1, newNeighbours2);
		if (compare!=0){
			return compare;
		}

    	List<List<List<Atom>>> prioritisedNeighbours1 = formListsWithSamePriority(newNeighbours1);
    	List<List<List<Atom>>> prioritisedNeighbours2 = formListsWithSamePriority(newNeighbours2);

    	for (int i = 1; i <= prioritisedNeighbours1.size(); i++) {
    		List<List<Atom>> nextNeighbourLists1 = prioritisedNeighbours1.get(prioritisedNeighbours1.size() -i);
    		List<List<Atom>> nextNeighbourLists2 = prioritisedNeighbours2.get(prioritisedNeighbours2.size() -i);
	    	for (int j = 1; j <= nextNeighbourLists1.size(); j++) {
	    		List<Atom> nextNeighbours1 = nextNeighbourLists1.get(nextNeighbourLists1.size() -j);
	    		List<Atom> nextNeighbours2 = nextNeighbourLists2.get(nextNeighbourLists2.size() -j);
	    		List<Atom> newPreviousAtoms1 = new ArrayList<Atom>();
	    		for (Atom atom : nextNeighbours1) {
	    			newPreviousAtoms1.add(currentToPrevious1.get(atom));
				}
	    		List<Atom> newPreviousAtoms2 = new ArrayList<Atom>();
	    		for (Atom atom : nextNeighbours2) {
	    			newPreviousAtoms2.add(currentToPrevious2.get(atom));
				}
	    		CipState newCIPstate = new CipState(cipState.previousAtomToVisitedAtoms1, cipState.previousAtomToVisitedAtoms2, newPreviousAtoms1, newPreviousAtoms2, nextNeighbours1, nextNeighbours2);
	    		queue.add(newCIPstate);
	    	}
		}
    	return 0;
	}
	
	private int compareNeighboursByCIPpriorityRules(List<List<List<Atom>>> neighbours1, List<List<List<Atom>>> neighbours2){
    	int neighbours1Size = neighbours1.size();
    	int neighbours2Size = neighbours2.size();
    	int differenceInSize = neighbours1Size - neighbours2Size;
    	int maxCommonSize = neighbours1Size > neighbours2Size ? neighbours2Size : neighbours1Size;
    	for (int i = 1; i <= maxCommonSize; i++) {
			int difference = listOfAtomListsCIPComparator.compare(neighbours1.get(neighbours1Size -i), neighbours2.get(neighbours2Size -i));
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

	private List<List<List<Atom>>> getNextLevelNeighbours(Map<Atom,Set<Atom>> previousAtomToVisitedAtoms, List<Atom> previousAtoms, List<Atom> nextAtoms, Map<Atom, Atom> currentToPrevious ) {
		List<List<List<Atom>>> neighbours = getNextAtomsConvertingRevisitedToGhosts(nextAtoms, previousAtomToVisitedAtoms, previousAtoms, currentToPrevious);
		for (List<List<Atom>> list : neighbours) {
			Collections.sort(list, atomListCIPComparator);
		}
		Collections.sort(neighbours, listOfAtomListsCIPComparator);
		return neighbours;
	}
	
	
	/**
	 * If given say [H,C,C] this becomes [H] [C,C] 
	 * If given say [H,C,C] [H,C,C] this becomes [H,H] [C,C,C,C]
	 * If given say [H,C,C] [H,C,F] this becomes [H],[C,C][H][C][F]
	 * as [H,C,F] is higher priority than [H,C,C] so all its atoms must be evaluated first
	 * The original neighbours list is assumed to have been presorted.
	 * @param neighbours
	 */
	private List<List<List<Atom>>> formListsWithSamePriority(List<List<List<Atom>>> neighbours) {
		List<List<List<Atom>>> updatedNeighbours  = new LinkedList<List<List<Atom>>>();
		List<List<Atom>> listsToRemove  = new ArrayList<List<Atom>>();
		for (List<List<Atom>> neighbourLists : neighbours) {
			List<List<Atom>> updatedNeighbourLists  = new LinkedList<List<Atom>>();
			for (int i = 0; i < neighbourLists.size(); i++) {
				List<List<Atom>> neighbourListsToCombine = new ArrayList<List<Atom>>();
				List<Atom> primaryAtomList = neighbourLists.get(i);
				for (int j = i +1; j < neighbourLists.size(); j++) {
					if (atomListCIPComparator.compare(neighbourLists.get(i), neighbourLists.get(j))==0){
						neighbourListsToCombine.add(neighbourLists.get(j));
					}
				}
				for (List<Atom> neighbourList: neighbourListsToCombine) {
					listsToRemove.add(neighbourList);
					primaryAtomList.addAll(neighbourList);
				}
			}
			for (List<Atom> list : listsToRemove) {
				neighbourLists.remove(list);
			}
			//lists of same priority have been combined e.g. [H,C,C] [H,C,C] -->[H,C,C,H,C,C]
			for (int i = neighbourLists.size()-1; i >=0; i--) {
				List<Atom> neighbourList = neighbourLists.get(i);
				Collections.sort(neighbourList, atomicNumberComparator);
				int lastAtomicNumber = Integer.MAX_VALUE;
				List<Atom> currentAtomList = new ArrayList<Atom>();
				for (int j = neighbourList.size() -1; j >=0; j--) {
					Atom a = neighbourList.get(j);
					int atomicNumber = AtomProperties.elementToAtomicNumber.get(a.getElement());
					if (atomicNumber < lastAtomicNumber){
						if (!currentAtomList.isEmpty()){
							updatedNeighbourLists.add(0, currentAtomList);
						}
						currentAtomList =new ArrayList<Atom>();
						currentAtomList.add(a);
					}
					else{
						currentAtomList.add(a);
					}
					lastAtomicNumber = atomicNumber;
				}
				if (!currentAtomList.isEmpty()){
					updatedNeighbourLists.add(0, currentAtomList);
				}
			}
			updatedNeighbours.add(updatedNeighbourLists);
		}
		return updatedNeighbours;
	}


	/**
	 * Sorts atomLists by CIP rules, low to high
	 * @author dl387
	 *
	 */
	private class AtomListCIPComparator implements Comparator<List<Atom>> {
		public int compare(List<Atom> a, List<Atom> b){
	    	int aSize = a.size();
	    	int bSize = b.size();
	    	int differenceInSize = aSize - bSize;
	    	int maxCommonSize = aSize > bSize ? bSize : aSize;
	    	for (int i = 1; i <= maxCommonSize; i++) {
				int difference = compareByAtomicNumber(a.get(aSize -i), b.get(bSize -i));
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
	 * Sorts lists of atomLists by CIP rules, low to high
	 * @author dl387
	 *
	 */
	private class ListOfAtomListsCIPComparator implements Comparator<List<List<Atom>>> {
		public int compare(List<List<Atom>> a, List<List<Atom>> b){
	    	int aSize = a.size();
	    	int bSize = b.size();
	    	int differenceInSize = aSize - bSize;
	    	int maxCommonSize = aSize > bSize ? bSize : aSize;
	    	for (int i = 1; i <= maxCommonSize; i++) {
	    		List<Atom> aprime = a.get(aSize -i);
	    		List<Atom> bprime = b.get(bSize -i);
		    	int aprimeSize = aprime.size();
		    	int bprimeSize = bprime.size();
		    	int differenceInSizeprime = aprimeSize - bprimeSize;
		    	int maxCommonSizeprime = aprimeSize > bprimeSize ? bprimeSize : aprimeSize;
		    	for (int j = 1; j <= maxCommonSizeprime; j++) {
		    		int difference = compareByAtomicNumber(aprime.get(aprimeSize -j), bprime.get(bprimeSize -j));
					if (difference >0){
						return 1;
					}
					if (difference < 0){
						return -1;
					}
		    	}
		    	if (differenceInSizeprime >0){
		    		return 1;
		    	}
		    	if (differenceInSizeprime <0){
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
	
	private List<List<List<Atom>>> getNextAtomsConvertingRevisitedToGhosts(List<Atom> atoms, Map<Atom, Set<Atom>> previousAtomToVisitedAtoms, List<Atom> previousAtoms, Map<Atom, Atom> currentToPrevious) {
		List<List<List<Atom>>> allNeighbours = new ArrayList<List<List<Atom>>>();
		int counter =0;
		Atom lastPreviousAtom = null;
		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = atoms.get(i);
			Atom previousAtom = previousAtoms.get(i);
			List<Atom> neighbours = atom.getAtomNeighbours();
			neighbours.remove(previousAtoms.get(i));			
			replaceRevisitedAtomsWithGhosts(neighbours, previousAtomToVisitedAtoms.get(previousAtom), atom);
			previousAtomToVisitedAtoms.put(atom, new HashSet<Atom>(previousAtomToVisitedAtoms.get(previousAtom)));
			previousAtomToVisitedAtoms.get(atom).add(atom);
			Collections.sort(neighbours, atomicNumberComparator);
			if (lastPreviousAtom==null){
				lastPreviousAtom = previousAtom;
			}
			else if (lastPreviousAtom !=previousAtom){
				lastPreviousAtom = previousAtom;
				counter++;
			}
			if (allNeighbours.size() <=  counter){
				allNeighbours.add(new ArrayList<List<Atom>>());
			}
			allNeighbours.get(counter).add(neighbours);
			for (Atom neighbour : neighbours) {
				currentToPrevious.put(neighbour, atom);
			}
		}
		return allNeighbours;
	}

	private void replaceRevisitedAtomsWithGhosts(List<Atom> atoms, Set<Atom> visitedAtoms, Atom previousAtom) {
		for (int i = atoms.size() -1; i >=0; i--) {
			Atom neighbour = atoms.get(i);
			if (visitedAtoms.contains(neighbour)){//cycle detected
				atoms.remove(i);
				if (neighbour.getID()>0){//make sure not to ghost ghosts
					atoms.add(i, createGhostAtomFromAtom(neighbour, previousAtom));
				}
			}
		}
	}

	private Atom createGhostAtomFromAtom(Atom atomToGhost, Atom atomToBondItTo) {
		try{
			Atom ghost = new Atom(ghostIdCounter--, atomToGhost.getElement(), molecule);
			Bond b1 = new Bond(ghost, atomToBondItTo, 1);
			atomToBondItTo.addBond(b1);
			ghost.addBond(b1);
			molecule.addAtom(ghost);
			return ghost;
		} catch (StructureBuildingException e) {//should never throw unless molecule is null
			return null;
		}
	}
	
	/**
	 * Greater than 0 means a is preferred over b (vice versa for less than 1)
	 * @param a
	 * @param b
	 * @return
	 */
    private int compareByAtomicNumber(Atom a, Atom b){
    	int atomicNumber1 = AtomProperties.elementToAtomicNumber.get(a.getElement());
    	int atomicNumber2 = AtomProperties.elementToAtomicNumber.get(b.getElement());
    	if (atomicNumber1 > atomicNumber2){
    		return 1;
    	}
    	else if (atomicNumber1 < atomicNumber2){
    		return -1;
    	}
		return 0;
    }
    
	/**
	 * Removes the ghost atoms added by the CIP evaluating code
	 * @param idValue
	 */
	private void removeAnyGhostAtomsAddedAndCorrectGhostIdCounter(int idValue) {
		List<Atom> atomList = molecule.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getID() <= idValue){
				Atom adjacentAtom = atom.getAtomNeighbours().get(0);
				adjacentAtom.removeBond(atom.getFirstBond());
				molecule.removeAtom(atom);
			}
		}
		ghostIdCounter = idValue;
	}
}
