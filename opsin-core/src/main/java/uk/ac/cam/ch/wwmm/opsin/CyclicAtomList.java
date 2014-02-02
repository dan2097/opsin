package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

/**
 * Convenience class for iterating over a list of atoms that form a ring
 * Doing getNext when the index is the final atom in the list will return the first atom
 * Doing getPrevious when the index is the first atom in the list will return the final atom
 * @author dl387
 *
 */
class CyclicAtomList{
	private int index = -1;
	private final List<Atom> atomList;

	/**
	 * Construct a cyclicAtomList from an atomList
	 * Index defaults to -1
	 * @param atomList
	 */
	CyclicAtomList(List<Atom> atomList) {
		this.atomList = atomList;
	}
	
	/**
	 * Construct a cyclicAtomList from an atomList
	 * The second parameter sets the current index
	 * @param atomList
	 * @param index
	 * @throws StructureBuildingException
	 */
	CyclicAtomList(List<Atom> atomList, int index) throws StructureBuildingException {
		this.atomList = atomList;
		setIndex(index);
	}
	
	/**
	 * Returns the number of elements in this list. If this list contains more
	 * than <tt>Integer.MAX_VALUE</tt> elements, returns
	 * <tt>Integer.MAX_VALUE</tt>.
	 * 
	 * @return the number of elements in this list
	 */
	int size() {
		return atomList.size();
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
	 * Return the current index in the list
	 * @return
	 */
	int getIndex() {
		return index;
	}

	/**
	 * Set the current index
	 * @param index
	 * @throws StructureBuildingException
	 */
	void setIndex(int index) throws StructureBuildingException{
		if (index >= atomList.size()){
			throw new StructureBuildingException("Specified index is not within ringAtom list");
		}
		this.index = index;
	}

	/**
	 * Increments and returns the atom at the new index in the list (next atom)
	 * When the index is the final atom in the list will return the first atom
	 * @return
	 */
	Atom next() {
		int tempIndex = index + 1;
		if (tempIndex >= atomList.size()){
			tempIndex = 0;
		}
		index = tempIndex;
		return atomList.get(index);
	}
	
	/**
	 * Decrements and returns the atom at the new index in the list (previous atom)
	 * when the index is the first atom in the list will return the final atom
	 * @return
	 */
	Atom previous() {
		int tempIndex = index - 1;
		if (tempIndex < 0){
			tempIndex = atomList.size() -1 ;
		}
		index = tempIndex;
		return atomList.get(index);
	}
	
	/**
	 * Returns the next atom in the list
	 * When the index is the final atom in the list will return the first atom
	 * Doesn't effect the list
	 * @return
	 */
	Atom peekNext() {
		int tempIndex = index + 1;
		if (tempIndex >= atomList.size()){
			tempIndex = 0;
		}
		return atomList.get(tempIndex);
	}
	
	/**
	 * Returns the previous atom in the list
	 * when the index is the first atom in the list will return the final atom
	 * Doesn't effect the list
	 * @return
	 */
	Atom peekPrevious() {
		int tempIndex = index - 1;
		if (tempIndex < 0){
			tempIndex = atomList.size() -1 ;
		}
		return atomList.get(tempIndex);
	}

	/**
	 * Returns the atom corresponding to the current index
	 * Note that CycliAtomLists have a default index of -1
	 * @return
	 */
	Atom getCurrent() {
		return atomList.get(index);
	}
}
