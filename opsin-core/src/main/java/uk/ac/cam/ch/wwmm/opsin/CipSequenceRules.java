package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.List;
import java.util.Queue;

/**
 * An implementation of Rule 1 of the CIP rules i.e. constitutional differences excluding isotopes
 * Cases that require rules 2-5 to distinguish result in an exception
 * 
 * Phantom atoms are not added as I believe that the results of the program will still be the same even in their absence as everything beats a phantom and comparing phantoms to phantoms achieves nothing
 * (higher ligancy beats lower ligancy when comparisons are performed)
 * @author dl387
 *
 */
class CipSequenceRules {
	private static class CipOrderingRunTimeException extends RuntimeException {
		private static final long serialVersionUID = 1L;
		public CipOrderingRunTimeException(String message) {
			super(message);
		}
	}
	
	private final Atom chiralAtom;
	
    CipSequenceRules(Atom chiralAtom) {
		this.chiralAtom = chiralAtom;
	}
    
	/**
	 * Returns the chiral atom's neighbours in CIP order from lowest priority to highest priority
	 * @return
	 * @throws CipOrderingException 
	 */
	List<Atom> getNeighbouringAtomsInCIPOrder() throws CipOrderingException {
		List<Atom> neighbours = chiralAtom.getAtomNeighbours();
		try {
			Collections.sort(neighbours, new SortByCIPOrder(chiralAtom));
		}
		catch (CipOrderingRunTimeException e){
			throw new CipOrderingException(e.getMessage());
		}
		return neighbours;
	}
	
	/**
	 * Returns  the chiral atom's neighbours, with the exception of the given atom, in CIP order from lowest priority to highest priority
	 * @param neighbourToIgnore
	 * @return
	 * @throws CipOrderingException 
	 */
	List<Atom> getNeighbouringAtomsInCIPOrderIgnoringGivenNeighbour(Atom neighbourToIgnore) throws CipOrderingException {
		List<Atom> neighbours = chiralAtom.getAtomNeighbours();
		if (!neighbours.remove(neighbourToIgnore)){
			throw new IllegalArgumentException("OPSIN bug: Atom" + neighbourToIgnore.getID() +" was not a neighbour of the given stereogenic atom");
		}
		try {
			Collections.sort(neighbours, new SortByCIPOrder(chiralAtom));
		}
		catch (CipOrderingRunTimeException e){
			throw new CipOrderingException(e.getMessage());
		}
		return neighbours;
	}

	
	/**
	 * Holds information about what atoms to try next next and how those atoms were reached (to prevent immediate back tracking and to detect cycles)
	 * @author dl387
	 *
	 */
	private static class CipState{
		CipState(List<AtomWithHistory> nextAtoms1, List<AtomWithHistory> nextAtoms2) {
			this.nextAtoms1 = nextAtoms1;
			this.nextAtoms2 = nextAtoms2;
		}
		final List<AtomWithHistory> nextAtoms1;
		final List<AtomWithHistory> nextAtoms2;
	}
	
	/**
	 * Holds an atom with associated visited atoms
	 * @author dl387
	 *
	 */
	private static class AtomWithHistory{
		AtomWithHistory(Atom atom, List<Atom> visitedAtoms, Integer indexOfOriginalFromRoot) {
			this.atom = atom;
			this.visitedAtoms = visitedAtoms;
			this.indexOfOriginalFromRoot = indexOfOriginalFromRoot;
		}
		final Atom atom;
		final List<Atom> visitedAtoms;
		final Integer indexOfOriginalFromRoot;
	}
	
	/**
	 * Sorts atoms by their CIP order, low to high
	 * @author dl387
	 *
	 */
	private class SortByCIPOrder implements Comparator<Atom> {
		private final Atom chiralAtom;
		private final AtomListCIPComparator atomListCIPComparator = new AtomListCIPComparator();
		private final ListOfAtomListsCIPComparator listOfAtomListsCIPComparator = new ListOfAtomListsCIPComparator();
		private final CIPComparator cipComparator = new CIPComparator();
		private int rule = 0;
		

		SortByCIPOrder(Atom chiralAtom) {
			this.chiralAtom = chiralAtom;
		}
		
		public int compare(Atom a, Atom b){
	    	/*
	    	 * rule = 0 --> Rule 1a Higher atomic number precedes lower
	    	 * rule = 1 --> Rule 1b A duplicated atom, with its predecessor node having the same label closer to the root, ranks higher than a duplicated atom, with its predecessor node having the same label farther from the root, which ranks higher than any non-duplicated atom node
	    	 * rule = 2 --> Rule 2 Higher atomic mass number precedes lower
	    	 */
	    	for (rule = 0; rule <= 2; rule++) {
				List<Atom> atomsVisted = new ArrayList<Atom>();
				atomsVisted.add(chiralAtom);
				AtomWithHistory aWithHistory = new AtomWithHistory(a, atomsVisted, null);
				AtomWithHistory bWithHistory = new AtomWithHistory(b, new ArrayList<Atom>(atomsVisted), null);
				
	    		int compare = compareByCipRules(aWithHistory, bWithHistory);
				if (compare != 0){
					return compare;
				}
				
				List<AtomWithHistory> nextAtoms1 = new ArrayList<AtomWithHistory>();
				nextAtoms1.add(aWithHistory);
				
				List<AtomWithHistory> nextAtoms2 = new ArrayList<AtomWithHistory>();
				nextAtoms2.add(bWithHistory);

				CipState startingState = new CipState(nextAtoms1, nextAtoms2);
		    	Deque<CipState> cipStateQueue = new ArrayDeque<CipState>();
		    	cipStateQueue.add(startingState);
		    	/* Go through CIP states in a breadth-first manner:
		    	 * Neighbours of the given atom/s (if multiple atoms this is because so far the two paths leading to them have been equivalent) are evaluated for both a and b
		    	 * Neighbours are sorted by CIP priority
		    	 * Comparisons performed between neighbours of a and neighbours of b (will break if compare!=0)
		    	 * Degenerate neighbours grouped together
		    	 * CIP state formed for each list of neighbours and added to queue in order of priority
		    	 *
		    	 */
		    	while(!cipStateQueue.isEmpty()){
		    		CipState currentState = cipStateQueue.removeFirst();
		    		compare = compareAtNextLevel(currentState, cipStateQueue);
		    		if (compare != 0){
		    			return compare;
		    		}
		    	}
			}
	    	throw new CipOrderingRunTimeException("Failed to assign CIP stereochemistry, this indicates a bug in OPSIN or a limitation in OPSIN's implementation of the sequence rules");
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
			List<List<List<AtomWithHistory>>> newNeighbours1 = getNextLevelNeighbours(cipState.nextAtoms1);
			List<List<List<AtomWithHistory>>> newNeighbours2 = getNextLevelNeighbours(cipState.nextAtoms2);

			int compare = compareNeighboursByCIPpriorityRules(newNeighbours1, newNeighbours2);

			if (compare != 0){
				return compare;
			}
	    	List<List<List<AtomWithHistory>>> prioritisedNeighbours1 = formListsWithSamePriority(newNeighbours1);
	    	List<List<List<AtomWithHistory>>> prioritisedNeighbours2 = formListsWithSamePriority(newNeighbours2);

	    	//As earlier compare was 0,  prioritisedNeighbours/nextNeighbourLists1 and 2 are the same size
	    	for (int i = prioritisedNeighbours1.size() - 1; i >=0; i--) {
	     		List<List<AtomWithHistory>> nextNeighbourLists1 = prioritisedNeighbours1.get(i);
	    		List<List<AtomWithHistory>> nextNeighbourLists2 = prioritisedNeighbours2.get(i);
	    		for (int j = nextNeighbourLists1.size() - 1; j >=0; j--) {
		    		List<AtomWithHistory> nextNeighbours1 = nextNeighbourLists1.get(j);
		    		List<AtomWithHistory> nextNeighbours2 = nextNeighbourLists2.get(j);
		    		CipState newCIPstate = new CipState(nextNeighbours1, nextNeighbours2);
		    		queue.add(newCIPstate);
	    		}
			}
	    	return 0;
		}
		
		private int compareNeighboursByCIPpriorityRules(List<List<List<AtomWithHistory>>> neighbours1, List<List<List<AtomWithHistory>>> neighbours2){
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

		private List<List<List<AtomWithHistory>>> getNextLevelNeighbours(List<AtomWithHistory> nextAtoms) {
			List<List<List<AtomWithHistory>>> neighbours = getNextAtomsWithAppropriateGhostAtoms(nextAtoms);
			for (List<List<AtomWithHistory>> list : neighbours) {
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
		private List<List<List<AtomWithHistory>>> formListsWithSamePriority(List<List<List<AtomWithHistory>>> neighbours) {
			List<List<List<AtomWithHistory>>> updatedNeighbours  = new ArrayList<List<List<AtomWithHistory>>>();
			List<List<AtomWithHistory>> listsToRemove  = new ArrayList<List<AtomWithHistory>>();
			for (List<List<AtomWithHistory>> neighbourLists : neighbours) {
				for (int i = 0; i < neighbourLists.size(); i++) {
					List<List<AtomWithHistory>> neighbourListsToCombine = new ArrayList<List<AtomWithHistory>>();
					List<AtomWithHistory> primaryAtomList = neighbourLists.get(i);
					for (int j = i + 1; j < neighbourLists.size(); j++) {
						List<AtomWithHistory> neighbourListToCompareWith = neighbourLists.get(j);
						if (atomListCIPComparator.compare(primaryAtomList, neighbourListToCompareWith)==0){
							neighbourListsToCombine.add(neighbourListToCompareWith);
						}
					}
					for (List<AtomWithHistory> neighbourList: neighbourListsToCombine) {
						listsToRemove.add(neighbourList);
						primaryAtomList.addAll(neighbourList);
					}
				}
				neighbourLists.removeAll(listsToRemove);

				List<List<AtomWithHistory>> updatedNeighbourLists  = new ArrayList<List<AtomWithHistory>>();
				//lists of same priority have been combined e.g. [H,C,C] [H,C,C] -->[H,C,C,H,C,C]
				for (int i = 0, lnls = neighbourLists.size(); i < lnls; i++) {
					List<AtomWithHistory> neighbourList = neighbourLists.get(i);
					Collections.sort(neighbourList, cipComparator);
					AtomWithHistory lastAtom = null;
					List<AtomWithHistory> currentAtomList = new ArrayList<AtomWithHistory>();
					for (int j = 0, lnl = neighbourList.size(); j < lnl; j++) {
						AtomWithHistory a = neighbourList.get(j);
						if (lastAtom != null && compareByCipRules(lastAtom, a) != 0){
							updatedNeighbourLists.add(currentAtomList);
							currentAtomList = new ArrayList<AtomWithHistory>();
						}
						currentAtomList.add(a);
						lastAtom = a;
					}
					if (!currentAtomList.isEmpty()){
						updatedNeighbourLists.add(currentAtomList);
					}
				}
				updatedNeighbours.add(updatedNeighbourLists);
			}
			return updatedNeighbours;
		}


		/**
		 * Sorts atoms by their atomic number, low to high
		 * @author dl387
		 *
		 */
		private class CIPComparator implements Comparator<AtomWithHistory> {
		    public int compare(AtomWithHistory a, AtomWithHistory b){
		    	return compareByCipRules(a, b);
		    }
		}

		/**
		 * Sorts atomLists by CIP rules, low to high
		 * @author dl387
		 *
		 */
		private class AtomListCIPComparator implements Comparator<List<AtomWithHistory>> {
			public int compare(List<AtomWithHistory> a, List<AtomWithHistory> b){
		    	int aSize = a.size();
		    	int bSize = b.size();
		    	int differenceInSize = aSize - bSize;
		    	int maxCommonSize = aSize > bSize ? bSize : aSize;
		    	for (int i = 1; i <= maxCommonSize; i++) {
					int difference = compareByCipRules(a.get(aSize -i), b.get(bSize -i));
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
		private class ListOfAtomListsCIPComparator implements Comparator<List<List<AtomWithHistory>>> {
			public int compare(List<List<AtomWithHistory>> a, List<List<AtomWithHistory>> b){
		    	int aSize = a.size();
		    	int bSize = b.size();
		    	int differenceInSize = aSize - bSize;
		    	int maxCommonSize = aSize > bSize ? bSize : aSize;
		    	for (int i = 1; i <= maxCommonSize; i++) {
		    		List<AtomWithHistory> aprime = a.get(aSize -i);
		    		List<AtomWithHistory> bprime = b.get(bSize -i);
			    	int aprimeSize = aprime.size();
			    	int bprimeSize = bprime.size();
			    	int differenceInSizeprime = aprimeSize - bprimeSize;
			    	int maxCommonSizeprime = aprimeSize > bprimeSize ? bprimeSize : aprimeSize;
			    	for (int j = 1; j <= maxCommonSizeprime; j++) {
			    		int difference = compareByCipRules(aprime.get(aprimeSize -j), bprime.get(bprimeSize -j));
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
		
		/**
		 * Gets the neighbouring atoms bar the previous atoms
		 * If the neighbouring atom has already been visited it is replaced with a ghost atom
		 * Multiple bonds including those to previous atoms yield ghost atoms unless the bond goes to the chiral atom e.g. in a sulfoxide
		 * @param atoms
		 * @return
		 */
		private List<List<List<AtomWithHistory>>> getNextAtomsWithAppropriateGhostAtoms(List<AtomWithHistory> atoms) {
			List<List<List<AtomWithHistory>>> allNeighbours = new ArrayList<List<List<AtomWithHistory>>>();
			int counter =0;
			Atom lastPreviousAtom = null;
			for (AtomWithHistory atomWithHistory : atoms) {
				Atom atom = atomWithHistory.atom;
				List<Atom> visitedAtoms = atomWithHistory.visitedAtoms;
				Atom previousAtom = visitedAtoms.get(visitedAtoms.size()-1);
				List<Atom> visitedAtomsIncludingCurrentAtom = new ArrayList<Atom>(visitedAtoms);
				visitedAtomsIncludingCurrentAtom.add(atom);

				List<AtomWithHistory> neighboursWithHistory = new ArrayList<AtomWithHistory>();
				for(Bond b :  atom.getBonds()) {
					Atom atomBondConnectsTo = b.getOtherAtom(atom);
					if (!atomBondConnectsTo.equals(chiralAtom)){//P-91.1.4.2.4 (higher order bonds to chiral centre do not involve duplication of atoms)
						for (int j = b.getOrder(); j >1; j--) {//add ghost atoms to represent higher order bonds
							Atom ghost = new Atom(atomBondConnectsTo.getElement());
							if (rule > 0){
								int indexOfOriginalAtom = visitedAtoms.indexOf(atomBondConnectsTo);
								if (indexOfOriginalAtom != -1){
									neighboursWithHistory.add(new AtomWithHistory(ghost, visitedAtomsIncludingCurrentAtom, indexOfOriginalAtom));
								}
								else{
									neighboursWithHistory.add(new AtomWithHistory(ghost, visitedAtomsIncludingCurrentAtom, visitedAtoms.size() + 1));
								}
							}
							else{
								neighboursWithHistory.add(new AtomWithHistory(ghost, visitedAtomsIncludingCurrentAtom, null));
							}
						}
					}
					if (!atomBondConnectsTo.equals(previousAtom)){
						if (visitedAtoms.contains(atomBondConnectsTo)){//cycle detected, add ghost atom instead
							Atom ghost = new Atom(atomBondConnectsTo.getElement());
							if (rule > 0){
								neighboursWithHistory.add(new AtomWithHistory(ghost, visitedAtomsIncludingCurrentAtom, visitedAtoms.indexOf(atomBondConnectsTo)));
							}
							else{
								neighboursWithHistory.add(new AtomWithHistory(ghost, visitedAtomsIncludingCurrentAtom, null));
							}
						}
						else{
							neighboursWithHistory.add(new AtomWithHistory(atomBondConnectsTo, visitedAtomsIncludingCurrentAtom, null));
						}
					}
				}
				Collections.sort(neighboursWithHistory, cipComparator);
				if (lastPreviousAtom==null){
					lastPreviousAtom = previousAtom;
				}
				else if (lastPreviousAtom !=previousAtom){
					lastPreviousAtom = previousAtom;
					counter++;
				}
				if (allNeighbours.size() <=  counter){
					allNeighbours.add(new ArrayList<List<AtomWithHistory>>());
				}
				allNeighbours.get(counter).add(neighboursWithHistory);
			}
			return allNeighbours;
		}
		
		/**
		 * Greater than 0 means a is preferred over b (vice versa for less than 1)
		 * @param a
		 * @param b
		 * @return
		 */
	    private int compareByCipRules(AtomWithHistory a, AtomWithHistory b){
	    	//rule 1a
	    	//prefer higher atomic number
	    	int atomicNumber1 = a.atom.getElement().ATOMIC_NUM;
	    	int atomicNumber2 = b.atom.getElement().ATOMIC_NUM;
	    	if (atomicNumber1 > atomicNumber2){
	    		return 1;
	    	}
	    	else if (atomicNumber1 < atomicNumber2){
	    		return -1;
	    	}
	    	if (rule > 0){
	    		//rule 1b
	    		//prefer duplicate to non-duplicate
	    		Integer indexFromRoot1 = a.indexOfOriginalFromRoot;
	       		Integer indexFromRoot2 = b.indexOfOriginalFromRoot;
	    		if (indexFromRoot1 != null && indexFromRoot2 == null){
	    			return 1;
	    		}
	    		if (indexFromRoot1 == null && indexFromRoot2 != null){
	    			return -1;
	    		}
	    		//prefer duplicate of node closer to root
	    		if (indexFromRoot1 != null && indexFromRoot2 != null){
	    	 		if (indexFromRoot1 < indexFromRoot2 ){
		    			return 1;
		    		}
	    	 		if (indexFromRoot1 > indexFromRoot2 ){
		    			return -1;
		    		}
	    		}
	    		if (rule > 1){
		    		//rule 2
		    		//prefer higher atomic mass
	    	    	Integer atomicMass1 = a.atom.getIsotope();
	    	    	Integer atomicMass2 = b.atom.getIsotope();
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
	    		}
	    		
	    	}
			return 0;
	    }
	}

}
