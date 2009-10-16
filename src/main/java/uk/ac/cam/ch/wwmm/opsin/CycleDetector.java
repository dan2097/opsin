package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

class CycleDetector {
	
	/**
	 * Performs a depth first search for rings hence assigning whether atoms are in rings or not
	 * This is necessary for deciding the applicability, and in some cases meaning, of suffixes and to determine what atoms are capable of having spare valency
	 * Fragments made of disconnected sections are supported
	 * @param frag
	 * @throws StructureBuildingException
	 */
	static void assignWhetherAtomsAreInCycles(Fragment frag) throws StructureBuildingException{
		List<Atom> atomList =frag.getAtomList();
		HashSet<Atom> atomsVisited = new HashSet<Atom>();
		for (Atom a : atomList) {//as OPSIN does not disallow disconnected sections within a single "fragment" (e.g. in suffixes) for vigorousness this for loop is required
			if (atomsVisited.contains(a)){continue;}//typically for all but the first atom this will be true
			LinkedList<Atom> orderAtomsVisited = new LinkedList<Atom>();
			LinkedList<Atom[]> atomStack =new LinkedList<Atom[]>();

			//Start from the first atom in the atomList and add neighbours to the relevant lists
			List<Atom> neighbours =a.getAtomNeighbours();
			atomsVisited.add(a);
			orderAtomsVisited.add(a);
			for (Atom neighbour : neighbours) {
				atomStack.add(new Atom[]{neighbour,a});//record both the atom and the atom from which you came from
			}
	
			while (atomStack.size()>0){
				Atom[] nextAtomArray =atomStack.removeLast();//depth first traversal
				a =nextAtomArray[0];
				if (atomsVisited.contains(a)){continue;}//parent atom must have been in a ring and this atom already visited another way
				atomsVisited.add(a);
				orderAtomsVisited.add(a);
				neighbours =a.getAtomNeighbours();
				int stackSize =atomStack.size();
				for (Atom neighbour : neighbours) {
					if (neighbour.equals(nextAtomArray[1])){continue;}//the neighbour is the parent of the atom e.g. in CCC you went you went from the first C then tried to go to the first as it is a neighbour
					if (atomsVisited.contains(neighbour)){//atom has already been visited, this must be because a ring has been formed
						int i =orderAtomsVisited.indexOf(neighbour);
						for (; i < orderAtomsVisited.size(); i++) {//go through atoms in orderAtomsVisited from the previously visited neighbour up to the end of the list and mark all these atoms as being in a ring
							Atom ringAtom = orderAtomsVisited.get(i);
							ringAtom.setAtomIsInACycle(true);
						}
					}
					else{//add atom to stack, its neighbours will be recursively investigated shortly
						atomStack.add(new Atom[]{neighbour, a});
					}
				}
				if (stackSize==atomStack.size() && stackSize > 0){//stack has not been changed, so you must be at a terminus. If the stack is empty this step is not necessary
					Atom atomToBackTrackTo =atomStack.getLast()[1];
					while (!atomToBackTrackTo.equals(orderAtomsVisited.getLast())){//remove atoms until you get back to the point from which you are investigating atoms. I believe atoms removed in this way are acyclic atoms
						orderAtomsVisited.removeLast();
					}
				}
			}
		}
	}
}
