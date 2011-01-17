package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

class CycleDetector {

	/**
	 * Performs a depth first search for rings hence assigning whether atoms are in rings or not
	 * This is necessary for deciding the applicability, and in some cases meaning, of suffixes and to determine what atoms are capable of having spare valency
	 * Fragments made of disconnected sections are supported
	 * @param frag
	 */
	static void assignWhetherAtomsAreInCycles(Fragment frag) {
		List<Atom> atomList =frag.getAtomList();
		for (Atom atom : atomList) {
			atom.setAtomIsInACycle(false);
			atom.setProperty(Atom.VISITED, null);
		}
		for (Atom a : atomList) {//as OPSIN does not disallow disconnected sections within a single "fragment" (e.g. in suffixes) for vigorousness this for loop is required
			if(a.getProperty(Atom.VISITED)==null){//typically for all but the first atom this will be true
				traverseRings(a, null, 0);
			}
		}
	}
	
	private static int traverseRings(Atom currentAtom, Atom previousAtom, int depth){
		int temp, result;
		if(currentAtom.getProperty(Atom.VISITED)!=null){
			return currentAtom.getProperty(Atom.VISITED);
		}
		currentAtom.setProperty(Atom.VISITED, depth);
		result = depth+1;
		List<Atom> neighbours = currentAtom.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			if (neighbour.equals(previousAtom)){
				continue;
			}
			temp = traverseRings(neighbour, currentAtom, depth+1);
			if( temp <= depth) {
				result = Math.min(result,temp);
			}
		}
		if( result <= depth ){
			currentAtom.setAtomIsInACycle(true);
		}
		return result;

	}

	private static class IntraFragmentPathSearchState{
		final Atom currentAtom;
		final LinkedList<Atom> orderAtomsVisited;
		public IntraFragmentPathSearchState(Atom currentAtom, LinkedList<Atom> orderAtomsVisited ) {
			this.currentAtom = currentAtom;
			this.orderAtomsVisited = orderAtomsVisited;
		}
		Atom getCurrentAtom() {
			return currentAtom;
		}
		LinkedList<Atom> getOrderAtomsVisited() {
			return orderAtomsVisited;
		}
	}
	
	static List<List<Atom>> getIntraFragmentPathsBetweenAtoms(Atom a1, Atom a2, Fragment fragment) throws StructureBuildingException {
		List<List<Atom>> paths = new ArrayList<List<Atom>>();
		LinkedList<IntraFragmentPathSearchState> stateStack = new LinkedList<IntraFragmentPathSearchState>();
		stateStack.add(new IntraFragmentPathSearchState(a1, new LinkedList<Atom>()));
		while (stateStack.size()>0){
			IntraFragmentPathSearchState state  =stateStack.removeLast();//depth first traversal
			LinkedList<Atom> orderAtomsVisited = state.getOrderAtomsVisited();
			Atom nextAtom = state.getCurrentAtom();
			orderAtomsVisited.add(nextAtom);
			List<Atom> neighbours = fragment.getIntraFragmentAtomNeighbours(nextAtom);
			for (Atom neighbour : neighbours) {
				if (orderAtomsVisited.contains(neighbour)){//atom already visited by this path
					continue;
				}
				if (neighbour ==a2 ){//target atom found
					paths.add(new ArrayList<Atom>(orderAtomsVisited.subList(1, orderAtomsVisited.size())));
				}
				else{//add atom to stack, its neighbours will be recursively investigated shortly
					stateStack.add(new IntraFragmentPathSearchState(neighbour, new LinkedList<Atom>(orderAtomsVisited)));
				}
			}
		}
		return paths;
	}
}
