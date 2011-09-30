package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

/**
 * Convenience class for iterating over a list of atoms that form a ring
 * Doing getNext when the indice is the final atom in the list will return the first atom
 * Doing getPrevious when the indice is the first atom in the list will return the final atom
 * @author dl387
 *
 */
class CyclicAtomList{
	private int indice = -1;
	private final List<Atom> atomList;

	/**
	 * Construct a cyclicAtomList from an atomList
	 * Indice defaults to -1
	 * @param atomList
	 */
	CyclicAtomList(List<Atom> atomList) {
		this.atomList = atomList;
	}
	
	/**
	 * Construct a cyclicAtomList from an atomList
	 * The second parameter sets the current indice
	 * @param atomList
	 * @param indice
	 * @throws StructureBuildingException
	 */
	CyclicAtomList(List<Atom> atomList, int indice) throws StructureBuildingException {
		this.atomList = atomList;
		setIndice(indice);
	}

	/**
	 * Returns the atom at the specified position in this list.
	 * @param index index of the element to return
	 * @return Atom the atom at the specified position in this list
	 * @throws IndexOutOfBoundsException - if the index is out of range (index < 0 || index >= size())
	 */
	Atom get(int index) throws IndexOutOfBoundsException {
		return atomList.get(index);
	}

	/**
	 * Return the current indice in the list
	 * @return
	 */
	int getIndice() {
		return indice;
	}

	/**
	 * Set the current indice
	 * @param indice
	 * @throws StructureBuildingException
	 */
	void setIndice(int indice) throws StructureBuildingException{
		if (indice >= atomList.size()){
			throw new StructureBuildingException("Specified indice is not within ringAtom list");
		}
		this.indice =indice;
	}

	/**
	 * Returns the next atom in the list
	 * When the indice is the final atom in the list will return the first atom
	 * @return
	 */
	Atom getNext() {
		int tempIndice = indice + 1;
		if (tempIndice >= atomList.size()){
			tempIndice=0;
		}
		indice =tempIndice;
		return atomList.get(indice);
	}
	
	/**
	 * Returns the previous atom in the list
	 * when the indice is the first atom in the list will return the final atom
	 * @return
	 */
	Atom getPrevious() {
		int tempIndice = indice - 1;
		if (tempIndice < 0){
			tempIndice = atomList.size() -1 ;
		}
		indice =tempIndice;
		return atomList.get(indice);
	}

	/**
	 * Returns the atom corresponding to the current indice
	 * Note that CycliAtomLists have a default indice of -1
	 * @return
	 */
	Atom getCurrent() {
		return atomList.get(indice);
	}
}
