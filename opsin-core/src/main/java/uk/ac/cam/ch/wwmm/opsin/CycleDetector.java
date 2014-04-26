package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Assigns whether atoms are in rings or not
 * @author dl387
 *
 */
class CycleDetector {

	/**
	 * Performs a depth first search for rings hence assigning whether atoms are in rings or not
	 * This is necessary for deciding the applicability, and in some cases meaning, of suffixes and to determine what atoms are capable of having spare valency
	 * Fragments made of disconnected sections are supported
	 * @param frag
	 */
	static void assignWhetherAtomsAreInCycles(Fragment frag) {
		List<Atom> atomList = frag.getAtomList();
		for (Atom atom : atomList) {
			atom.setAtomIsInACycle(false);
			atom.setProperty(Atom.VISITED, null);
		}
		for (Atom a : atomList) {//as OPSIN does not disallow disconnected sections within a single "fragment" (e.g. in suffixes) for vigorousness this for loop is required
			if(a.getProperty(Atom.VISITED) == null){//true for only the first atom in a fully connected molecule
				traverseRings(a, null, 0);
			}
		}
	}
	
	private static int traverseRings(Atom currentAtom, Atom previousAtom, int depth){
		Integer previouslyAssignedDepth = currentAtom.getProperty(Atom.VISITED);
		if(previouslyAssignedDepth != null){
			return previouslyAssignedDepth;
		}
		currentAtom.setProperty(Atom.VISITED, depth);
		List<Atom> equivalentAtoms = new ArrayList<Atom>();
		equivalentAtoms.add(currentAtom);
		
		List<Atom> neighbours;
		for(;;) {
			//Non-recursively process atoms in a chain
			//add the atoms in the chain to equivalentAtoms as either all or none of them are in a ring
			neighbours = currentAtom.getAtomNeighbours();
			neighbours.remove(previousAtom);
			if (neighbours.size() != 1) {
				break;
			}
			Atom nextAtom = neighbours.get(0);
			if (nextAtom.getProperty(Atom.VISITED) != null) {
				//chain reached a previously visited atom, must be a ring
				break;
			}
			previousAtom = currentAtom;
			currentAtom = nextAtom;
			equivalentAtoms.add(currentAtom);
			currentAtom.setProperty(Atom.VISITED, ++depth);
		}

		int result = depth + 1;
		for (Atom neighbour : neighbours) {
		  int temp = traverseRings(neighbour, currentAtom, depth + 1);
		  result = Math.min(result, temp);
		}
		if (result < depth){
			for (Atom a : equivalentAtoms) {
				a.setAtomIsInACycle(true);
			}
		} else if (result == depth) {
			currentAtom.setAtomIsInACycle(true);
		}
		return result;
	}

	private static class PathSearchState{
		final Atom currentAtom;
		final List<Atom> orderAtomsVisited;
		public PathSearchState(Atom currentAtom, List<Atom> orderAtomsVisited ) {
			this.currentAtom = currentAtom;
			this.orderAtomsVisited = orderAtomsVisited;
		}
		Atom getCurrentAtom() {
			return currentAtom;
		}
		List<Atom> getOrderAtomsVisited() {
			return orderAtomsVisited;
		}
	}
	
	/**
	 * Attempts to find paths from a1 to a2 using only the given bonds
	 * @param a1
	 * @param a2
	 * @param peripheryBonds
	 * @return
	 */
	static List<List<Atom>> getPathBetweenAtomsUsingBonds(Atom a1, Atom a2, Set<Bond> peripheryBonds){
		List<List<Atom>> paths = new ArrayList<List<Atom>>();
		Deque<PathSearchState> stateStack = new ArrayDeque<PathSearchState>();
		stateStack.add(new PathSearchState(a1, new ArrayList<Atom>()));
		while (stateStack.size()>0){
			PathSearchState state  =stateStack.removeLast();//depth first traversal
			List<Atom> orderAtomsVisited = state.getOrderAtomsVisited();
			Atom nextAtom = state.getCurrentAtom();
			orderAtomsVisited.add(nextAtom);
			Set<Bond> neighbourBonds = new LinkedHashSet<Bond>(nextAtom.getBonds());
			neighbourBonds.retainAll(peripheryBonds);
			for (Bond neighbourBond : neighbourBonds) {
				Atom neighbour = neighbourBond.getOtherAtom(nextAtom);
				if (orderAtomsVisited.contains(neighbour)){//atom already visited by this path
					continue;
				}
				if (neighbour ==a2 ){//target atom found
					paths.add(new ArrayList<Atom>(orderAtomsVisited.subList(1, orderAtomsVisited.size())));
				}
				else{//add atom to stack, its neighbours will be recursively investigated shortly
					stateStack.add(new PathSearchState(neighbour, new ArrayList<Atom>(orderAtomsVisited)));
				}
			}
		}
		return paths;
	}
}
