package uk.ac.cam.ch.wwmm.opsin;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import nu.xom.Element;


/**A "struct" to hold the results of fragment building.*/
class BuildResults {
	/**Holds the ID or IDs of atoms that are currently marked as radicals.
	 * Typically these will be utilised by a word rule e.g. the ethyl of ethyl ethanoate has one
	 * Also holds the order of the bond that will be created when it is used (valency)
	 * setExplicitly says whether the outID absolutely definitely refers to that id or not.
	 * e.g. propyl is stored as prop-1-yl with this set to false while prop-2-yl has it set to true
	 * These OutIDs are the same objects as are present in the fragments*/
	private LinkedList<OutID> outIDs;

	/**The ID or IDs of the atoms that may be used to from things like esters*/
	private LinkedList<FunctionalID> functionalIDs;

	/**The ID or IDs of the atoms that must have bonds formed to them. Rarely used expect in the root of multiplicative names*/
	private LinkedList<InID> inIDs;

	/**A list of fragments that have been evaluated to form this BuildResults. They are in the order they would be found in the XML*/
	private LinkedHashSet<Fragment> fragments;

	/**A BuildResults is constructed from a list of Fragments.
	 * This constructor creates this list from the groups present in an XML word/bracket/sub element.*/
	BuildResults(BuildState state, Element wordSubOrBracket) {
		outIDs = new LinkedList<OutID>();
		functionalIDs = new LinkedList<FunctionalID>();
		inIDs = new LinkedList<InID>();
		fragments = new LinkedHashSet<Fragment>();
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(wordSubOrBracket, "group");
		for (Element group : groups) {
			Fragment frag = state.xmlFragmentMap.get(group);
			fragments.add(frag);
			outIDs.addAll(frag.getOutIDs());
			functionalIDs.addAll(frag.getFunctionalIDs());
			inIDs.addAll(frag.getInIDs());
		}
	}

	/**
	 * Construct a blank buildResults
	 */
	BuildResults() {
		outIDs = new LinkedList<OutID>();
		functionalIDs = new LinkedList<FunctionalID>();
		inIDs = new LinkedList<InID>();
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
	 * Returns the atom corresponding to position i in the outID list
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException
	 */
	Atom getOutAtom(int i) throws StructureBuildingException {
		OutID outID = outIDs.get(i);
		return outID.frag.getAtomByIDOrThrow(outID.id);
	}

	/**
	 * Returns the atom corresponding to position i in the outIDs list
	 * If not set explicitly and atom would violate valency or break aromaticity another is looked for
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException
	 */
	Atom getOutAtomTakingIntoAccountWhetherSetExplicitly(int i) throws StructureBuildingException {
		OutID outID = outIDs.get(i);
		if (outID.setExplicitly){
			return getOutAtom(i);
		}
		else{
			return outID.frag.getAtomByIdOrNextSuitableAtomOrThrow(outID.id, outID.valency);
		}
	}

	OutID getOutID(int i){
		return outIDs.get(i);
	}

	int getOutIDCount(){
		return outIDs.size();
	}

	void removeOutID(int i) throws StructureBuildingException{
		OutID outID =outIDs.get(i);
		if (outID.frag!=null){
			outID.frag.removeOutID(outID);
		}
		outIDs.remove(i);
	}

	void removeAllOutIDs() throws StructureBuildingException{
		for (int i = outIDs.size() -1; i >=0 ; i--) {
			removeOutID(i);
		}
	}

	/**
	 * Returns the atom corresponding to position i in the functionalIDs list
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException
	 */
	Atom getFunctionalAtom(int i) throws StructureBuildingException {
		FunctionalID functionalId =functionalIDs.get(i);
		return functionalId.frag.getAtomByIDOrThrow(functionalId.id);
	}

	void removeFunctionalID(int i) {
		FunctionalID functionalID =functionalIDs.get(i);
		if (functionalID.frag!=null){
			functionalID.frag.removeFunctionalID(functionalID);
		}
		functionalIDs.remove(i);
	}

	int getFunctionalIDCount(){
		return functionalIDs.size();
	}

	/**
	 * Returns the first OutId
	 * @return OutID
	 */
	OutID getFirstOutID(){
		return outIDs.get(0);
	}

	int getInIDCount(){
		return inIDs.size();
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
		outIDs.addAll(otherBR.outIDs);
		functionalIDs.addAll(otherBR.functionalIDs);
		inIDs.addAll(otherBR.inIDs);
		fragments.addAll(otherBR.fragments);
	}
}
