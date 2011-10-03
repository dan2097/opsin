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

public class CipSequenceRules {
	
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
	    	if (compare==0){
	    		throw new RuntimeException("Failed to assign CIP stereochemistry, this indicates a bug in OPSIN or a limitation in OPSIN's implementation of the sequence rules");
	    	}
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
			
			Set<Atom> atomsVisted = new HashSet<Atom>();
			atomsVisted.add(chiralAtom);
			Map<Atom,Set<Atom>> previousAtomToVisitedAtoms1 = new HashMap<Atom, Set<Atom>>();
			previousAtomToVisitedAtoms1.put(chiralAtom, atomsVisted);
			Map<Atom,Set<Atom>> previousAtomToVisitedAtoms2 = new HashMap<Atom, Set<Atom>>();
			previousAtomToVisitedAtoms2.put(chiralAtom, new HashSet<Atom>(atomsVisted));

	    	CipState startingState = new CipState(previousAtomToVisitedAtoms1, previousAtomToVisitedAtoms2, previousAtoms1, previousAtoms2, nextAtoms1, nextAtoms2);
			return startingState;
		}
	}
	
	private final static AtomListCIPComparator atomListCIPComparator = new AtomListCIPComparator();
	private final static ListOfAtomListsCIPComparator listOfAtomListsCIPComparator = new ListOfAtomListsCIPComparator();
	private final static AtomicNumberComparator atomicNumberComparator = new AtomicNumberComparator();
	private final Fragment molecule;
	private int ghostIdCounter =-1;
	//phantom atoms are not added as I believe that the results of the program will still be the same even in their absence as everything beats a phantom and comparing phantoms to phantoms achieves nothing

    CipSequenceRules(Fragment molecule) {
		this.molecule = molecule;
	}
    
	/**
	 * Returns the given atoms neighbours in CIP order from lowest priority to highest priority
	 * @param chiralAtom
	 * @return
	 */
	List<Atom> getNeighbouringAtomsInCIPOrder(Atom chiralAtom) {
		List<Atom> neighbours = chiralAtom.getAtomNeighbours();
		addGhostAtomsForCIPAssignment(chiralAtom);
		Collections.sort(neighbours, new SortByCIPOrder(chiralAtom));
		removeGhostAtoms();
		return neighbours;
	}
	
	/**
	 * Returns the given atoms neighbours, with the exception of the given atom, in CIP order from lowest priority to highest priority
	 * @param chiralAtom
	 * @return
	 */
	List<Atom> getNeighbouringAtomsInCIPOrderIgnoringGivenNeighbour(Atom chiralAtom, Atom neighbourToIgnore) {
		List<Atom> neighbours = chiralAtom.getAtomNeighbours();
		if (!neighbours.remove(neighbourToIgnore)){
			throw new RuntimeException("OPSIN bug: " + neighbourToIgnore.toCMLAtom().toXML() +" was not a neighbour of the given stereogenic atom");
		}
		addGhostAtomsForCIPAssignment(chiralAtom);
		Collections.sort(neighbours, new SortByCIPOrder(chiralAtom));
		removeGhostAtoms();
		return neighbours;
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
	private static class AtomListCIPComparator implements Comparator<List<Atom>> {
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
	private static class ListOfAtomListsCIPComparator implements Comparator<List<List<Atom>>> {
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
		Atom ghost = new Atom(ghostIdCounter--, atomToGhost.getElement(), molecule);
		Bond b1 = new Bond(ghost, atomToBondItTo, 1);
		atomToBondItTo.addBond(b1);
		ghost.addBond(b1);
		molecule.addAtom(ghost);
		return ghost;
	}
	
	/**
	 * Greater than 0 means a is preferred over b (vice versa for less than 1)
	 * @param a
	 * @param b
	 * @return
	 */
    private static int compareByAtomicNumber(Atom a, Atom b){
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
	
	/**
	 * Adds "ghost" atoms in accordance with the CIP rules for handling double bonds
	 * e.g. C=C --> C(G)=C(G) where ghost is a carbon with no hydrogen bonded to it
	 * Higher order bonds connected to the chiral atom are not converted in accordance with P-91.1.4.2.4 (IUPAC 2004 guidelines)
	 * @param chiralAtom Higher order bonds connected to this atom are not touched
	 */
	private void addGhostAtomsForCIPAssignment(Atom chiralAtom) {
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
	 * Removes the ghost atoms added by addGhostAtomsForCIPAssignment
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
}
