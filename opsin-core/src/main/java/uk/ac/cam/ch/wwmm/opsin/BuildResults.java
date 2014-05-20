package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * A "struct" to hold the results of fragment building.
 * @author dl387
 *
 */
class BuildResults {
	/**Holds the atoms that are currently marked as radicals. An atom may be listed twice for say diyl
	 * Typically these will be utilised by a word rule e.g. the ethyl of ethyl ethanoate has one
	 * Also holds the order of the bond that will be created when it is used (valency)
	 * setExplicitly says whether the outAtom absolutely definitely refers to that atom or not.
	 * e.g. propyl is stored as prop-1-yl with this set to false while prop-2-yl has it set to true
	 * These OutAtoms are the same objects as are present in the fragments*/
	private final List<OutAtom> outAtoms;

	/**The atoms that may be used to from things like esters*/
	private final List<FunctionalAtom> functionalAtoms;

	/**A list of fragments that have been evaluated to form this BuildResults. They are in the order they would be found in the XML*/
	private final Set<Fragment> fragments;

	/**A BuildResults is constructed from a list of Fragments.
	 * This constructor creates this list from the groups present in an XML word/bracket/sub element.
     * @param wordSubOrBracket*/
	BuildResults(Element wordSubOrBracket) {
		outAtoms = new ArrayList<OutAtom>();
		functionalAtoms = new ArrayList<FunctionalAtom>();
		fragments = new LinkedHashSet<Fragment>();
		List<Element> groups = OpsinTools.getDescendantElementsWithTagName(wordSubOrBracket, XmlDeclarations.GROUP_EL);
		for (Element group : groups) {
			Fragment frag = group.getFrag();
			fragments.add(frag);
			for (int i = 0, l = frag.getOutAtomCount(); i < l; i++) {
				outAtoms.add(frag.getOutAtom(i));
			}
			int functionalAtomCount = frag.getFunctionalAtomCount();
			if (functionalAtomCount > 0){
				Element parent = group.getParent();
				if (parent.getName().equals(XmlDeclarations.ROOT_EL) ||
						OpsinTools.getNextGroup(group) == null){
					for (int i = 0, l = functionalAtomCount; i < l; i++) {
						functionalAtoms.add(frag.getFunctionalAtom(i));
					}
				}
			}
			
		}
	}

	/**
	 * Construct a blank buildResults
	 */
	BuildResults() {
		outAtoms = new ArrayList<OutAtom>();
		functionalAtoms = new ArrayList<FunctionalAtom>();
		fragments = new LinkedHashSet<Fragment>();
	}

	/**
	 * Returns a read only view of the fragments in this BuildResults
	 * @return
	 */
	Set<Fragment> getFragments(){
		return Collections.unmodifiableSet(fragments);
	}

	int getFragmentCount(){
		return fragments.size();
	}

	/**
	 * Returns the atom corresponding to position i in the outAtoms list
	 * If not set explicitly and atom would violate valency or break aromaticity another is looked for
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException
	 */
	Atom getOutAtomTakingIntoAccountWhetherSetExplicitly(int i) throws StructureBuildingException {
		OutAtom outAtom = outAtoms.get(i);
		if (outAtom.isSetExplicitly()){
			return outAtom.getAtom();
		}
		else{
			return outAtom.getAtom().getFrag().getAtomOrNextSuitableAtomOrThrow(outAtom.getAtom(), outAtom.getValency(), false);
		}
	}

	OutAtom getOutAtom(int i){
		return outAtoms.get(i);
	}

	int getOutAtomCount(){
		return outAtoms.size();
	}

	OutAtom removeOutAtom(int i) {
		OutAtom outAtom =outAtoms.get(i);
		outAtom.getAtom().getFrag().removeOutAtom(outAtom);
		return outAtoms.remove(i);
	}

	void removeAllOutAtoms() {
		for (int i = outAtoms.size() -1; i >=0 ; i--) {
			removeOutAtom(i);
		}
	}

	/**
	 * Returns the atom corresponding to position i in the functionalAtoms list
	 * @param i index
	 * @return atom
	 */
	Atom getFunctionalAtom(int i){
		return functionalAtoms.get(i).getAtom();
	}

	FunctionalAtom removeFunctionalAtom(int i) {
		FunctionalAtom functionalAtom =functionalAtoms.get(i);
		functionalAtom.getAtom().getFrag().removeFunctionalAtom(functionalAtom);
		return functionalAtoms.remove(i);
	}

	int getFunctionalAtomCount(){
		return functionalAtoms.size();
	}

	/**
	 * Returns the first OutAtom
	 * @return OutAtom
	 */
	OutAtom getFirstOutAtom(){
		return outAtoms.get(0);
	}

	/**
	 * Returns the atom corresponding to the given id assuming the atom the id corresponds to is within the list of fragment in this Buildresults
	 * @param id index
	 * @return atom
	 * @throws StructureBuildingException
	 */
	Atom getAtomByIdOrThrow(int id) throws StructureBuildingException {
		for (Fragment fragment : fragments) {
			Atom outAtom =fragment.getAtomByID(id);
			if (outAtom!=null){
				return outAtom;
			}
		}
		throw new StructureBuildingException("No fragment contained this id: " + id);
	}

	void mergeBuildResults(BuildResults otherBR) {
		outAtoms.addAll(otherBR.outAtoms);
		functionalAtoms.addAll(otherBR.functionalAtoms);
		fragments.addAll(otherBR.fragments);
	}

	/**
	 * Returns the sum of the charges of the fragments in the buildResults
	 * @return
	 */
	int getCharge() {
		int totalCharge=0;
		for (Fragment frag : fragments) {
			totalCharge+=frag.getCharge();
		}
		return totalCharge;
	}
}
