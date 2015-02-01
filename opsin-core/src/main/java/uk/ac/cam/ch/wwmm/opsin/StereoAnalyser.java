package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

/**
 * Identifies stereocentres and determines the CIP order of connected atoms
 * @author dl387
 *
 */
class StereoAnalyser {
	/** The atoms/bonds upon which this StereoAnalyser is operating */
	private final Collection<Atom> atoms;
	private final Collection<Bond> bonds;
	
	/** Maps each atom to its currently assigned colour. Eventually all atoms in non identical environments will have different colours. Higher is higher priority*/
	private final Map<Atom, Integer> mappingToColour;

	/** Maps each atom to an array of the colours of its neighbours*/
	private final Map<Atom, int[]> atomNeighbourColours;

	private final AtomNeighbouringColoursComparator atomNeighbouringColoursComparator = new AtomNeighbouringColoursComparator();
	private final static AtomicNumberThenAtomicMassComparator atomicNumberThenAtomicMassComparator = new AtomicNumberThenAtomicMassComparator();
	
	/**
	 * Holds information about a tetrahedral stereocentre
	 * @author dl387
	 *
	 */
	class StereoCentre{
		private final Atom stereoAtom;
		private final boolean trueStereoCentre;

		/**
		 * Creates a stereocentre object from a tetrahedral stereocentre atom
		 * @param stereoAtom
		 */
		StereoCentre(Atom stereoAtom, Boolean isTrueStereoCentre) {
			this.stereoAtom = stereoAtom;
			this.trueStereoCentre = isTrueStereoCentre;
		}

		Atom getStereoAtom() {
			return stereoAtom;
		}
		
		/**
		 * Does this atom have 4 constitutionally different groups (or 3 and a lone pair)
		 * or is it only a stereo centre due to the presence of other centres in the molecule
		 * @return
		 */
		boolean isTrueStereoCentre() {
			return trueStereoCentre;
		}

		List<Atom> getCipOrderedAtoms() throws CipOrderingException {
			List<Atom> cipOrderedAtoms = new CipSequenceRules(stereoAtom).getNeighbouringAtomsInCIPOrder();
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
		 * @throws CipOrderingException 
		 */
		List<Atom> getOrderedStereoAtoms() throws CipOrderingException {
			Atom a1 = bond.getFromAtom();
			Atom a2 = bond.getToAtom();
			List<Atom> cipOrdered1 = new CipSequenceRules(a1).getNeighbouringAtomsInCIPOrderIgnoringGivenNeighbour(a2);
			List<Atom> cipOrdered2 = new CipSequenceRules(a2).getNeighbouringAtomsInCIPOrderIgnoringGivenNeighbour(a1);
			List<Atom> stereoAtoms = new ArrayList<Atom>();
			stereoAtoms.add(cipOrdered1.get(cipOrdered1.size()-1));//highest CIP adjacent to a1
			stereoAtoms.add(a1);
			stereoAtoms.add(a2);
			stereoAtoms.add(cipOrdered2.get(cipOrdered2.size()-1));//highest CIP adjacent to a2
			return stereoAtoms;
		}
	}
	
	/**
	 * Sorts atoms by their atomic number, low to high
	 * In the case of a tie sorts by atomic mass
	 * @author dl387
	 *
	 */
	private static class AtomicNumberThenAtomicMassComparator implements Comparator<Atom> {
	    public int compare(Atom a, Atom b){
	    	return compareAtomicNumberThenAtomicMass(a, b);
	    }
	}
	
	private static int compareAtomicNumberThenAtomicMass(Atom a, Atom b){
    	int atomicNumber1 = a.getElement().ATOMIC_NUM;
    	int atomicNumber2 = b.getElement().ATOMIC_NUM;
    	if (atomicNumber1 > atomicNumber2){
    		return 1;
    	}
    	else if (atomicNumber1 < atomicNumber2){
    		return -1;
    	}
    	Integer atomicMass1 = a.getIsotope();
    	Integer atomicMass2 = b.getIsotope();
    	if (atomicMass1 != null && atomicMass2 == null){
    		return 1;
    	}
    	else if (atomicMass1 == null && atomicMass2 != null){
    		return -1;
    	}
    	else if (atomicMass1 != null && atomicMass2 != null){
        	if (atomicMass1 > atomicMass2){
	    		return 1;
	    	}
	    	else if (atomicMass1 < atomicMass2){
	    		return -1;
	    	}
    	}
		return 0;
	}
	
	/**
	 * Sorts based on the list of colours for neighbouring atoms 
	 * e.g. [1,2] > [1,1]  [1,1,3] > [2,2,2]  [1,1,3] > [3]  
	 * @author dl387
	 *
	 */
	private class AtomNeighbouringColoursComparator implements Comparator<Atom> {
	    public int compare(Atom a, Atom b){
	    	int[] colours1 = atomNeighbourColours.get(a);
	    	int[] colours2 = atomNeighbourColours.get(b);
	    	
	    	int colours1Size = colours1.length;
	    	int colours2Size = colours2.length;
	    	
	    	int maxCommonColourSize = Math.min(colours1Size, colours2Size);
	    	for (int i = 1; i <= maxCommonColourSize; i++) {
				int difference = colours1[colours1Size - i] - colours2[colours2Size - i];
				if (difference > 0){
					return 1;
				}
				if (difference < 0){
					return -1;
				}
			}
	    	int differenceInSize = colours1Size - colours2Size;
	    	if (differenceInSize > 0){
	    		return 1;
	    	}
	    	if (differenceInSize < 0){
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
	 */
	StereoAnalyser(Fragment molecule) {
		this (molecule.getAtomList(), molecule.getBondSet());
	}

	/**
	 * Employs a derivative of the InChI algorithm to label which atoms are equivalent.
	 * These labels can then be used by the findStereo(Atoms/Bonds) functions to find features that
	 * can possess stereoChemistry
	 * NOTE: All bonds of every atom must be in the set of bonds, no atom may have a bond to an atom not in the list
	 * @param atoms
	 * @param bonds
	 */
	StereoAnalyser(Collection<Atom> atoms, Collection<Bond> bonds) {
		this.atoms = atoms;
		this.bonds = bonds;
		List<Atom> ghostAtoms = addGhostAtoms();
		List<Atom> atomsToSort = new ArrayList<Atom>(atoms);
		atomsToSort.addAll(ghostAtoms);
		mappingToColour = new HashMap<Atom, Integer>(atomsToSort.size());
		atomNeighbourColours = new HashMap<Atom, int[]>(atomsToSort.size());
		Collections.sort(atomsToSort, atomicNumberThenAtomicMassComparator);
		List<List<Atom>> groupsByColour = populateColoursByAtomicNumberAndMass(atomsToSort);
		boolean changeFound = true;
		while(changeFound){
			for (List<Atom> groupWithAColour : groupsByColour) {
				for (Atom atom : groupWithAColour) {
					int[] neighbourColours = findColourOfNeighbours(atom);
					atomNeighbourColours.put(atom, neighbourColours);
				}
			}
			List<List<Atom>> updatedGroupsByColour = new ArrayList<List<Atom>>();
			changeFound = populateColoursAndReportIfColoursWereChanged(groupsByColour, updatedGroupsByColour);
			groupsByColour = updatedGroupsByColour;
		}
		removeGhostAtoms(ghostAtoms);
	}

	/**
	 * Adds "ghost" atoms in the same way as the CIP rules for handling double bonds
	 * e.g. C=C --> C(G)=C(G) where ghost is a carbon with no hydrogen bonded to it
	 * @return The ghost atoms created
	 */
	private List<Atom> addGhostAtoms() {
		List<Atom> ghostAtoms = new ArrayList<Atom>();
		for (Bond bond : bonds) {
			int bondOrder = bond.getOrder();
			for (int i = bondOrder; i > 1; i--) {
				Atom fromAtom = bond.getFromAtom();
				Atom toAtom = bond.getToAtom();

				Atom ghost1 = new Atom(fromAtom.getElement());
				Bond b1 = new Bond(ghost1, toAtom, 1);
				toAtom.addBond(b1);
				ghost1.addBond(b1);
				ghostAtoms.add(ghost1);
				
				Atom ghost2 = new Atom(toAtom.getElement());
				Bond b2 = new Bond(ghost2, fromAtom, 1);
				fromAtom.addBond(b2);
				ghost2.addBond(b2);
				ghostAtoms.add(ghost2);
			}
		}
		return ghostAtoms;
	}

	/**
	 * Removes the ghost atoms added by addGhostAtoms
	 * @param ghostAtoms 
	 */
	private void removeGhostAtoms(List<Atom> ghostAtoms) {
		for (Atom atom : ghostAtoms) {
			Bond b = atom.getFirstBond();
			Atom adjacentAtom = b.getOtherAtom(atom);
			adjacentAtom.removeBond(atom.getFirstBond());
		}
	}


	/**
	 * Takes a list of atoms sorted by atomic number/mass
	 * and populates the mappingToColour map
	 * @param atomList
	 * @return 
	 */
	private List<List<Atom>> populateColoursByAtomicNumberAndMass(List<Atom> atomList) {
		List<List<Atom>> groupsByColour = new ArrayList<List<Atom>>();
		Atom previousAtom = null;
		List<Atom> atomsOfThisColour = new ArrayList<Atom>();
		int atomsSeen = 0;
		for (Atom atom : atomList) {
			if (previousAtom != null && compareAtomicNumberThenAtomicMass(previousAtom, atom) != 0){
				for (Atom atomOfthisColour : atomsOfThisColour) {
					mappingToColour.put(atomOfthisColour, atomsSeen);
				}
				groupsByColour.add(atomsOfThisColour);
				atomsOfThisColour = new ArrayList<Atom>();
			}
			previousAtom = atom;
			atomsOfThisColour.add(atom);
			atomsSeen++;
		}
		if (!atomsOfThisColour.isEmpty()){
			for (Atom atomOfThisColour : atomsOfThisColour) {
				mappingToColour.put(atomOfThisColour, atomsSeen);
			}
			groupsByColour.add(atomsOfThisColour);
		}
		return groupsByColour;
	}

	/**
	 * Takes the lists of atoms pre-grouped by colour and sorts each by its neighbours colours
	 * The updatedGroupsByColour is populated with those for which this process caused a change
	 * and populates the mappingToColour map
	 * Returns whether mappingToColour was changed 
	 * @param groupsByColour 
	 * @param updatedGroupsByColour 
	 * @return boolean Whether mappingToColour was changed
	 */
	private boolean populateColoursAndReportIfColoursWereChanged(List<List<Atom>> groupsByColour, List<List<Atom>> updatedGroupsByColour) {
		boolean changeFound = false;
		int atomsSeen = 0;
		for (List<Atom> groupWithAColour : groupsByColour) {
			Collections.sort(groupWithAColour, atomNeighbouringColoursComparator);
			Atom previousAtom = null;
			List<Atom> atomsOfThisColour = new ArrayList<Atom>();
			for (Atom atom : groupWithAColour) {
				if (previousAtom != null && atomNeighbouringColoursComparator.compare(previousAtom, atom) != 0){
					for (Atom atomOfThisColour : atomsOfThisColour) {
						if (!changeFound && atomsSeen != mappingToColour.get(atomOfThisColour)){
							changeFound = true;
						}
						mappingToColour.put(atomOfThisColour, atomsSeen);
					}
					updatedGroupsByColour.add(atomsOfThisColour);
					atomsOfThisColour = new ArrayList<Atom>();
				}
				previousAtom = atom;
				atomsOfThisColour.add(atom);
				atomsSeen++;
			}
			if (!atomsOfThisColour.isEmpty()){
				for (Atom atomOfThisColour : atomsOfThisColour) {
					if (!changeFound && atomsSeen != mappingToColour.get(atomOfThisColour)){
						changeFound = true;
					}
					mappingToColour.put(atomOfThisColour, atomsSeen);
				}
				updatedGroupsByColour.add(atomsOfThisColour);
			}
		}
		return changeFound;
	}

	/**
	 * Produces a sorted (low to high) array of the colour of the atoms surrounding a given atom
	 * @param atom
	 * @return int[] colourOfAdjacentAtoms
	 */
	private int[] findColourOfNeighbours(Atom atom) {	
		List<Bond> bonds = atom.getBonds();
		int bondCount = bonds.size();
		int[] colourOfAdjacentAtoms = new int[bondCount];
		for (int i = 0; i < bondCount; i++) {
			Bond bond = bonds.get(i);
			Atom otherAtom = bond.getOtherAtom(atom);
			colourOfAdjacentAtoms[i] = mappingToColour.get(otherAtom);
		} 
		Arrays.sort(colourOfAdjacentAtoms);//sort such that this goes from low to high
		return colourOfAdjacentAtoms;
	}

	/**
	 * Retrieves a list of any tetrahedral stereoCentres
	 * Internally this is done by checking whether the "colour" of all neighbouring atoms of the tetrahedral atom are different
	 * @return List<StereoCentre>
	 */
	List<StereoCentre> findStereoCentres(){
		List<Atom> potentialStereoAtoms = getPotentialStereoCentres();
		List<Atom> trueStereoCentres = new ArrayList<Atom>();
		for (Atom potentialStereoAtom : potentialStereoAtoms) {
			if (isTrueStereCentre(potentialStereoAtom)){
				trueStereoCentres.add(potentialStereoAtom);
			}
		}
		List<StereoCentre> stereoCentres = new ArrayList<StereoCentre>();
		for (Atom trueStereoCentreAtom : trueStereoCentres) {
			stereoCentres.add(new StereoCentre(trueStereoCentreAtom, true));
		}

		potentialStereoAtoms.removeAll(trueStereoCentres);
		List<Atom> paraStereoCentres = findParaStereoCentres(potentialStereoAtoms, trueStereoCentres);
		for (Atom paraStereoCentreAtom : paraStereoCentres) {
			stereoCentres.add(new StereoCentre(paraStereoCentreAtom, false));
		}
		return stereoCentres;
	}

	/**
	 * Retrieves atoms that pass the isPossiblyStereogenic() criteria
	 * @return
	 */
	private List<Atom> getPotentialStereoCentres() {
		List<Atom> potentialStereoAtoms = new ArrayList<Atom>();
		for (Atom atom : atoms) {
			if (isPossiblyStereogenic(atom)){
				potentialStereoAtoms.add(atom);
			}
		}
		return potentialStereoAtoms;
	}
	
	/**
	 * Checks whether the atom has 3 or 4 neighbours all of which are constitutionally different
	 * @param potentialStereoAtom
	 * @return
	 */
	private boolean isTrueStereCentre(Atom potentialStereoAtom) {
		List<Atom> neighbours = potentialStereoAtom.getAtomNeighbours();
		if (neighbours.size()!=3 && neighbours.size()!=4){
			return false;
		}
		int[] colours = new int[4];
		for (int i = neighbours.size() - 1 ; i >=0; i--) {
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
		return !foundIdenticalNeighbour;
	}
	
	/**
	 * Finds a subset of the stereocentres associated with rule 2 from:
	 * DOI: 10.1021/ci00016a003
	 * @param potentialStereoAtoms
	 * @param trueStereoCentres 
	 */
	private List<Atom> findParaStereoCentres(List<Atom> potentialStereoAtoms, List<Atom> trueStereoCentres) {
		List<Atom> paraStereoCentres = new ArrayList<Atom>();
		for (Atom potentialStereoAtom : potentialStereoAtoms) {
			List<Atom> neighbours = potentialStereoAtom.getAtomNeighbours();
			if (neighbours.size() == 4){
				int[] colours = new int[4];
				for (int i = neighbours.size() - 1 ; i >=0; i--) {
					colours[i] = mappingToColour.get(neighbours.get(i));
				}
				//find pairs of constitutionally identical substituents
				Map<Integer, Integer> foundPairs = new HashMap<Integer, Integer>();
				for (int i = 0; i < 4; i++) {
					int cl = colours[i];
					for (int j = i +1; j < 4; j++) {
						if (cl == colours[j]){
							foundPairs.put(i, j);
							break;
						}
					}
				}
				int pairs = foundPairs.keySet().size();
				if (pairs==1 || pairs==2){
					boolean hasTrueStereoCentreInAllBranches = true;
					for (Entry<Integer, Integer> entry: foundPairs.entrySet()) {
						if (!branchesHaveTrueStereocentre(neighbours.get(entry.getKey()), neighbours.get(entry.getValue()), potentialStereoAtom, trueStereoCentres)){
							hasTrueStereoCentreInAllBranches = false;
							break;
						}
					}
					if (hasTrueStereoCentreInAllBranches){
						paraStereoCentres.add(potentialStereoAtom);
					}
				}
			}
		}
		return paraStereoCentres;
	}


	private boolean branchesHaveTrueStereocentre(Atom branchAtom1, Atom branchAtom2, Atom potentialStereoAtom, List<Atom> trueStereoCentres) {
		List<Atom> atomsToVisit= new ArrayList<Atom>();
		Set<Atom> visitedAtoms = new HashSet<Atom>();
		visitedAtoms.add(potentialStereoAtom);
		atomsToVisit.add(branchAtom1);
		atomsToVisit.add(branchAtom2);
		while(!atomsToVisit.isEmpty()){
			List<Atom> newAtomsToVisit = new ArrayList<Atom>();
			while(!atomsToVisit.isEmpty()){
				Atom atom = atomsToVisit.remove(0);
				if (trueStereoCentres.contains(atom)){
					return true;
				}
				if (atomsToVisit.contains(atom)){//the two branches have converged on this atom, don't investigate neighbours of it
					do{
						atomsToVisit.remove(atom);
					}
					while (atomsToVisit.contains(atom));
					continue;
				}
				else{
					List<Atom> neighbours = atom.getAtomNeighbours();
					for (Atom neighbour : neighbours) {
						if (visitedAtoms.contains(neighbour)){
							continue;
						}
						newAtomsToVisit.add(neighbour);
					}
				}
				visitedAtoms.add(atom);
			}
			atomsToVisit = newAtomsToVisit;
		}
		return false;
	}

	/**
	 * Checks whether an atom could be a tetrahedral stereocentre by checking that it is both tetrahedral
	 * and does not have neighbours that are identical due to resonance/tautomerism
	 * @param atom
	 * @return
	 */
	static boolean isPossiblyStereogenic(Atom atom){
		return isKnownPotentiallyStereogenic(atom) && !isAchiralDueToResonanceOrTautomerism(atom);
	}

	/**
	 * Roughly corresponds to the list of atoms in table 8 of the InChI manual
	 * Essentially does a crude check for whether an atom is known to be able to possess tetrahedral geometry
	 * and whether it is currently tetrahedral. Atoms that are tetrahedral but not typically considered chiral
	 * like tertiary amines are not recognised
	 * @param atom
	 * @return
	 */
	static boolean isKnownPotentiallyStereogenic(Atom atom) {
		List<Atom> neighbours = atom.getAtomNeighbours();
		ChemEl chemEl = atom.getElement();
		if (neighbours.size() == 4){
			if (chemEl == ChemEl.B || chemEl == ChemEl.C || chemEl == ChemEl.Si || chemEl == ChemEl.Ge ||
					chemEl == ChemEl.Sn || chemEl == ChemEl.N || chemEl == ChemEl.P || chemEl == ChemEl.As ||
						chemEl == ChemEl.S || chemEl == ChemEl.Se){
				return true;
			}
		}
		else if (neighbours.size() ==3){
			if ((chemEl == ChemEl.S || chemEl == ChemEl.Se) && (atom.getIncomingValency()==4 || (atom.getCharge() ==1 && atom.getIncomingValency()==3))){
				//tetrahedral sulfur/selenium - 3 bonds and the lone pair
				return true;
			}
			if (chemEl == ChemEl.N && atom.getCharge() ==0 && atom.getIncomingValency()==3 && atomsContainABondBetweenThemselves(neighbours)){
				return true;
				//nitrogen where two attached atoms are connected together
			}
		}
		return false;
	}
	
	private static boolean atomsContainABondBetweenThemselves(List<Atom> atoms) {
		for (Atom atom : atoms) {
			for (Atom neighbour : atom.getAtomNeighbours()) {
				if (atoms.contains(neighbour)){
					return true;
				}
			}
		}
		return false;
	}

	static boolean isAchiralDueToResonanceOrTautomerism(Atom atom) {
		ChemEl chemEl = atom.getElement();
		if(chemEl == ChemEl.N || 
				chemEl == ChemEl.P || 
				chemEl == ChemEl.As || 
				chemEl == ChemEl.S || 
				chemEl == ChemEl.Se) {
			List<Atom> neighbours = atom.getAtomNeighbours();
			Set<String> resonanceAndTautomerismAtomicElementPlusIsotopes = new HashSet<String>();
			for (Atom neighbour : neighbours) {
				ChemEl neighbourChemEl = neighbour.getElement();
				if ((neighbourChemEl.isChalcogen() || neighbourChemEl == ChemEl.N)
						&& isOnlyBondedToHydrogensOtherThanGivenAtom(neighbour, atom)){
					if (resonanceAndTautomerismAtomicElementPlusIsotopes.contains(neighbourChemEl.toString() + atom.getIsotope())){
						return true;
					}
					resonanceAndTautomerismAtomicElementPlusIsotopes.add(neighbourChemEl.toString() + atom.getIsotope());
				}
				if (neighbourChemEl == ChemEl.H && neighbour.getBondCount()==1){
					//terminal H atom neighbour
					return true;
				}
			}
		}
		return false;
	}

	private static boolean isOnlyBondedToHydrogensOtherThanGivenAtom(Atom atom, Atom attachedNonHydrogen) {
		for (Atom neighbour: atom.getAtomNeighbours()) {
			if (neighbour.equals(attachedNonHydrogen)){
				continue;
			}
			if (neighbour.getElement() != ChemEl.H){
				return false;
			}
		}
		return true;
	}

	/**
	 *  Retrieves a list of any double bonds possessing the potential to have E/Z stereoChemistry
	 *  This is done internally by checking the two atoms attached to the ends of the double bond are different
	 *  As an exception nitrogen's lone pair is treated like a low priority group and so is allowed to only have 1 atom connected to it
	 * @return
	 */
	List<StereoBond> findStereoBonds() {
		List<StereoBond> stereoBonds = new ArrayList<StereoBond>();
		for (Bond bond : bonds) {
			if (bond.getOrder()==2){
				Atom a1 = bond.getFromAtom();
				List<Atom> neighbours1 =  a1.getAtomNeighbours();
				neighbours1.remove(bond.getToAtom());
				if (neighbours1.size()==2 || (neighbours1.size()==1 && a1.getElement() == ChemEl.N && a1.getIncomingValency()==3 && a1.getCharge()==0)){
					if (neighbours1.size()==2 && mappingToColour.get(neighbours1.get(0)).equals(mappingToColour.get(neighbours1.get(1)))){
						continue;
					}
					Atom a2 = bond.getToAtom();
					List<Atom> neighbours2 = a2.getAtomNeighbours();
					neighbours2.remove(bond.getFromAtom());
					if (neighbours2.size()==2 || (neighbours2.size()==1 && a2.getElement() == ChemEl.N && a2.getIncomingValency()==3 && a2.getCharge()==0)){
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
	 * Returns a number describing the environment of an atom. Atoms with the same number are in identical environments
	 * Null if atom was not part of this environment analysis
	 * @param a
	 * @return
	 */
	Integer getAtomEnvironmentNumber(Atom a) {
		return mappingToColour.get(a);
	}
}
