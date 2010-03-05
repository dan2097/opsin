package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nu.xom.Element;


/** Holds the Fragments during the construction of the molecule,
 * handles the building of new fragments and handles the creation/deletion of atoms/bonds
 *
 * @author ptc24
 *
 */
class FragmentManager {
	/** All of the atom-containing fragments in the molecule */
	private final Set<Fragment> fragPile;
	/** All of the inter-fragment bonds */
	private final Set<Bond> bondPile;
	/** A mapping between fragments and inter fragment bonds */
	private final Map<Fragment,Set<Bond>> fragToInterFragmentBond;
	/** A builder for fragments specified as SMILES */
	private final SMILESFragmentBuilder sBuilder;
	/** A builder for fragments specified as references to a CML data file */
	private final CMLFragmentBuilder cmlBuilder;
	/** A source of unique integers */
	private final IDManager idManager;

	/** Sets up a new Fragment mananger, containing no fragments.
	 *
	 * @param sBuilder A SMILESFragmentBuilder - dependency injection.
	 * @param cmlBuilder A CMLFragmentBuilder - dependency injection.
	 * @param idManager An IDManager.
	 *
	 * @throws Exception If the CML fragment file can't be found or otherwise used
	 */
	FragmentManager(SMILESFragmentBuilder sBuilder, CMLFragmentBuilder cmlBuilder, IDManager idManager) {
		if (sBuilder == null || cmlBuilder == null || idManager == null ){
			throw new IllegalArgumentException("FragmentManager was parsed a null object in its constructor!");
		}
		this.sBuilder = sBuilder;
		this.cmlBuilder = cmlBuilder;
		this.idManager = idManager;
		fragPile = new LinkedHashSet<Fragment>();
		bondPile = new LinkedHashSet<Bond>();
		fragToInterFragmentBond = new HashMap<Fragment, Set<Bond>>();
	}

	/** Builds a fragment, based on a reference to a CML data file
	 *
	 * @param idStr The name of the fragment in the CML file
	 * @param type The fragment type
	 * @param subType The fragment subType
	 * @return The fragment
	 * @throws StructureBuildingException If the fragment can't be built
	 */
	Fragment buildCML(String idStr, String type, String subType) throws StructureBuildingException {
		Fragment newFrag = cmlBuilder.build(idStr, type, subType, this);
		addFragment(newFrag);
		return newFrag;
	}

	/** Builds a fragment, based on an SMILES string
	 *
	 * @param smiles The fragment to build
	 * @return The built fragment
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles) throws StructureBuildingException {
		return buildSMILES(smiles, "", "");
	}

	/** Builds a fragment, based on an SMILES string
	 *
	 * @param smiles The fragment to build
	 * @param type The fragment type
	 * @param labelMapping How to label the fragment
	 * @return The built fragment
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles, String type, String labelMapping) throws StructureBuildingException {
		return buildSMILES(smiles, type, "", labelMapping);
	}

	/** Builds a fragment, based on an SMILES string
	 *
	 * @param smiles The fragment to build
	 * @param type The fragment type
	 * @param subType The fragment subType
	 * @param labelMapping How to label the fragment
	 * @return The built fragment
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles, String type, String subType, String labelMapping) throws StructureBuildingException {
		Fragment newFrag = sBuilder.build(smiles, type, subType, labelMapping, this);
		addFragment(newFrag);
		return newFrag;
	}

	/**Creates a new fragment, containing all of the atoms and bonds
	 * of all of the other fragments - ie the whole molecule. This is non-destructive,
	 * and does not update which fragment the Atoms think they are in. Atoms and Bonds
	 * are not copied.
	 *
	 * @return The unified fragment
	 * @throws StructureBuildingException 
	 */
	Fragment getUnifiedFragment() throws StructureBuildingException {
		Fragment outFrag = new Fragment();
		addFragment(outFrag);
		List<Fragment> fragments = new ArrayList<Fragment>(fragPile);
		for(Fragment f : fragments) {
			incorporateFragment(f, outFrag);//merge all fragments into one
		}
		return outFrag;
	}

	/** Incorporates a fragment, usually a suffix, into a parent fragment
	 * This does:
	 * Imports all of the atoms and bonds from another fragment into this one.
	 * Also imports outIDs/inIDs and functionalIDs
	 * Reassigns inter fragment bonds of the parent fragment as either intra fragment bonds
	 * of the parent fragment or as inter fragment bonds of the parent fragment
	 *
	 * The original fragment still maintains its original atomList/bondList/interFragmentBondList which is necessary for stereochemistry handling
	 *
	 * @param childFrag The fragment to be incorporated
	 * @param parentFrag The parent fragment
	 */
	void incorporateFragment(Fragment childFrag, Fragment parentFrag) throws StructureBuildingException {
		for(Atom atom : childFrag.getAtomList()) {
			parentFrag.addAtom(atom);
		}
		for(Bond bond : childFrag.getBondSet()) {
			parentFrag.addBond(bond);
		}
		for (OutID outID: childFrag.getOutIDs()) {
			outID.frag =parentFrag;
		}
		for (InID inID: childFrag.getInIDs()) {
			inID.frag =parentFrag;
		}
		for (FunctionalID functionalID: childFrag.getFunctionalIDs()) {
			functionalID.frag =parentFrag;
		}
		parentFrag.addOutIDs(childFrag.getOutIDs());
		parentFrag.addInIDs(childFrag.getInIDs());
		parentFrag.addFunctionalIDs(childFrag.getFunctionalIDs());

		for (Bond bond : fragToInterFragmentBond.get(childFrag)) {//reassign inter fragment bonds of child
			if (bond.getFromAtom().getFrag() ==parentFrag || bond.getToAtom().getFrag() ==parentFrag){
				if (bond.getFromAtom().getFrag() ==parentFrag && bond.getToAtom().getFrag() ==parentFrag){
					//bond is now enclosed within parentFrag so make it an intra fragment bond
					//and remove it from the interfragment list of the parentFrag
					parentFrag.addBond(bond);
					fragToInterFragmentBond.get(parentFrag).remove(bond);
				}
				else{
					//bond was an interfragment bond between the childFrag and another frag
					//It is now between the parentFrag and another frag
					addInterFragmentBond(bond);
				}
			}
		}
		if (!fragPile.remove(childFrag)){
			throw new StructureBuildingException("Fragment not found in fragPile");
		}
	}
	
	/** Incorporates a fragment, usually a suffix, into a parent fragment, creating a bond between them.
	 *
	 * @param childFrag The fragment to be incorporated
	 * @param fromID An id on that fragment
	 * @param parentFrag The parent fragment
	 * @param toID An id on that fragment
	 * @param bondOrder The order of the joining bond
	 */
	void incorporateFragment(Fragment childFrag, int fromID, Fragment parentFrag, int toID, int bondOrder) throws StructureBuildingException {
		incorporateFragment(childFrag, parentFrag);
		createBond(parentFrag.getAtomByIDOrThrow(fromID), parentFrag.getAtomByIDOrThrow(toID), bondOrder);
	}

	/** Converts an atom in a fragment to a different atomic symbol.
	 * Charged atoms can also be specified using a SMILES formula eg. [N+]
	 *
	 * @param a The atom to change to a heteroatom
	 * @param atomSymbol The atomic symbol to be used
	 * @param assignLocant Whether a locant should be assigned to the heteroatom if the locant is not used elsewhere
	 * @throws StructureBuildingException if the atom could not be found
	 */
	void makeHeteroatom(Atom a, String atomSymbol, boolean assignLocant) throws StructureBuildingException {
		if(atomSymbol.startsWith("[")) {
			Fragment f = sBuilder.build(atomSymbol, this);
			Atom referenceAtom = f.getAtomList().get(0);
			atomSymbol =referenceAtom.getElement();
			a.setCharge(referenceAtom.getCharge());
		}
		a.setElement(atomSymbol);
		a.removeElementSymbolLocants();
		if (assignLocant && a.getFrag().getAtomByLocant(atomSymbol) ==null){//if none of that element currently present add element symbol locant
			a.addLocant(atomSymbol);
		}
	}

	/** Gets an atom, given an id number
	 * Use this if you don't know what fragment the atom is in
	 * @param id The id of the atom
	 * @return The atom, or null if no such atom exists.
	 */
	Atom getAtomByID(int id) {
		for(Fragment f : fragPile) {
			Atom a = f.getAtomByID(id);
			if(a != null) return a;
		}
		return null;
	}

	/** Gets an atom, given an id number, throwing if fails.
	 * Use this if you don't know what fragment the atom is in
	 * @param id The id of the atom
	 * @return The atom
	 */
	Atom getAtomByIDOrThrow(int id) throws StructureBuildingException {
		Atom a = getAtomByID(id);
		if(a == null) throw new StructureBuildingException("Couldn't get atom by id");
		return a;
	}

	/**Turns all of the spare valencies in the fragments into double bonds.
	 *
	 * @throws StructureBuildingException
	 */
	void convertSpareValenciesToDoubleBonds() throws StructureBuildingException {
		for(Fragment f : fragPile) {
			f.convertSpareValenciesToDoubleBonds();
		}
	}

	/**
	 * Checks valencies are all chemically reasonable. An exception is thrown if any are not
	 * @throws StructureBuildingException
	 */
	void checkValencies() throws StructureBuildingException {
		for(Fragment f : fragPile) {
			f.checkValencies();
		}
	}

	Set<Bond> getBondPile() {
		return Collections.unmodifiableSet(bondPile);
	}

	Set<Fragment> getFragPile() {
		return Collections.unmodifiableSet(fragPile);
	}

	/**
	 * Adds a fragment to the fragPile
	 * @param frag
	 */
	private void addFragment(Fragment frag)  {
		fragPile.add(frag);
		fragToInterFragmentBond.put(frag, new HashSet<Bond>());
	}

	/**
	 * Removes a fragment from the fragPile and inter fragment bonds associated with it from the bondpile/fragToInterFragmentBond.
	 * Throws an exception if fragment wasn't present
	 * @param frag
	 * @throws StructureBuildingException
	 */
	void removeFragment(Fragment frag) throws StructureBuildingException {
		if (!fragPile.remove(frag)){
			throw new StructureBuildingException("Fragment not found in fragPile");
		}
		List<Bond> interFragmentBondsInvolvingFragment = new ArrayList<Bond>(fragToInterFragmentBond.get(frag));
		for (Bond bond : interFragmentBondsInvolvingFragment) {
			if (bond.getFromAtom().getFrag() ==frag){
				fragToInterFragmentBond.get(bond.getToAtom().getFrag()).remove(bond);
			}
			else{
				fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).remove(bond);
			}
			bondPile.remove(bond);
		}
		fragToInterFragmentBond.get(frag).clear();
	}

	int getOverallCharge() {
		int totalCharge=0;
		for (Fragment frag : fragPile) {
			totalCharge+=frag.getCharge();
		}
		return totalCharge;
	}

	/**
	 * Creates a copy of a fragment by copying data
	 * labels the atoms using new ids from the idManager
	 * @param originalFragment
	 * @return the clone of the fragment
	 * @throws StructureBuildingException
	 */
	Fragment copyAndRelabel(Fragment originalFragment) throws StructureBuildingException {
		return copyAndRelabel(originalFragment, null);
	}


	/**
	 * Creates a copy of a fragment by copying data
	 * labels the atoms using new ids from the idManager
	 * @param originalFragment
	 * @param stringToAddToAllLocants: typically used to append primes to all locants, can be null
	 * @return the clone of the fragment
	 * @throws StructureBuildingException
	 */
	Fragment copyAndRelabel(Fragment originalFragment, String stringToAddToAllLocants) throws StructureBuildingException {
		Fragment newFragment =new Fragment(originalFragment.getType(), originalFragment.getSubType());
		HashMap<Integer, Integer> idMap = new HashMap<Integer, Integer>();//maps old ID to new ID
		List<Atom> atomList =originalFragment.getAtomList();
		newFragment.setIndicatedHydrogen(originalFragment.getIndicatedHydrogen());
		List<OutID> outIDs =originalFragment.getOutIDs();
		List<FunctionalID> functionalIDs =originalFragment.getFunctionalIDs();
		List<InID> inIDs =originalFragment.getInIDs();
		int defaultInId =originalFragment.getDefaultInID();
		for (Atom atom : atomList) {
			int id = idManager.getNextID();
			ArrayList<String> newLocants = new ArrayList<String>(atom.getLocants());
			if (stringToAddToAllLocants !=null){
				for (int i = 0; i < newLocants.size(); i++) {
					newLocants.set(i, newLocants.get(i) + stringToAddToAllLocants);
				}
			}
			Atom newAtom =new Atom(id, atom.getElement(), newFragment);
			for (String newLocant : newLocants) {
				newAtom.addLocant(newLocant);
			}
			newAtom.setCharge(atom.getCharge());
			newAtom.setSpareValency(atom.hasSpareValency());
			newAtom.setProtonsExplicitlyAddedOrRemoved(atom.getProtonsExplicitlyAddedOrRemoved());
			if (atom.getAtomParityElement() != null){
				newAtom.setAtomParityElement(new Element(atom.getAtomParityElement()));
			}
			newAtom.setExplicitHydrogens(atom.getExplicitHydrogens());
			newAtom.setLambdaConventionValency(atom.getLambdaConventionValency());
			//outValency is derived from the outIDs so is automatically cloned
			newAtom.setAtomIsInACycle(atom.getAtomIsInACycle());
			newAtom.setType(atom.getType());//may be different from fragment type if the original atom was formerly in a suffix
			newAtom.setNotes(new HashMap<String, String>(atom.getNotes()));
			newFragment.addAtom(newAtom);
			idMap.put(atom.getID(),id);
		}
        for (OutID outID : outIDs) {
            newFragment.addOutID(idMap.get(outID.id), outID.valency, outID.setExplicitly);
        }
		for (FunctionalID functionalID : functionalIDs) {
			newFragment.addFunctionalID(idMap.get(functionalID.id));
		}
        for (InID inID : inIDs) {
            newFragment.addInID(idMap.get(inID.id), inID.valency);
        }
		newFragment.setDefaultInID(idMap.get(defaultInId));
		Set<Bond> bondSet =originalFragment.getBondSet();
		for (Bond bond : bondSet) {
			Bond newBond = createBond(newFragment.getAtomByIDOrThrow(idMap.get(bond.getFrom())),newFragment.getAtomByIDOrThrow(idMap.get(bond.getTo())), bond.getOrder());
			newBond.setSmilesStereochemistry(bond.getSmilesStereochemistry());
			if (bond.getBondStereoElement() != null){
				newBond.setBondStereoElement(new Element(bond.getBondStereoElement()));
			}
		}
		addFragment(newFragment);
		return newFragment;
	}

	/**
	 * Takes an element and produces a copy of it. Groups and suffixes are copied so that the new element
	 * has its own group and suffix fragments
	 * @param elementToBeCloned
	 * @param state The current buildstate
	 * @return
	 * @throws StructureBuildingException
	 */
	Element cloneElement(BuildState state, Element elementToBeCloned) throws StructureBuildingException {
		return cloneElement(state, elementToBeCloned, "");
	}

	/**
	 * Takes an element and produces a copy of it. Groups and suffixes are copied so that the new element
	 * has its own group and suffix fragments
	 * @param elementToBeCloned
	 * @param state The current buildstate
	 * @param stringToAddToAllLocants A string to append to all locants in the cloned fragments
	 * @return
	 * @throws StructureBuildingException
	 */
	Element cloneElement(BuildState state, Element elementToBeCloned, String stringToAddToAllLocants) throws StructureBuildingException {
		Element clone = new Element(elementToBeCloned);
		List<Element> originalGroups = XOMTools.getDescendantElementsWithTagName(elementToBeCloned, "group");
		List<Element> clonedGroups = XOMTools.getDescendantElementsWithTagName(clone, "group");
		HashMap<Fragment,Fragment> oldNewFragmentMapping  =new HashMap<Fragment, Fragment>();
		for (int i = 0; i < originalGroups.size(); i++) {
			Fragment originalFragment =state.xmlFragmentMap.get(originalGroups.get(i));
			Fragment newFragment = copyAndRelabel(originalFragment, stringToAddToAllLocants);
			oldNewFragmentMapping.put(originalFragment, newFragment);
			state.xmlFragmentMap.put(clonedGroups.get(i), newFragment);
			ArrayList<Fragment> originalSuffixes =state.xmlSuffixMap.get(originalGroups.get(i));
			ArrayList<Fragment> newSuffixFragments =new ArrayList<Fragment>();
			for (Fragment suffix : originalSuffixes) {
				newSuffixFragments.add(state.fragManager.copyAndRelabel(suffix));
			}
			state.xmlSuffixMap.put(clonedGroups.get(i), newSuffixFragments);
		}
		Set<Bond> interFragmentBondsToClone = new HashSet<Bond>();
		for (Fragment originalFragment : oldNewFragmentMapping.keySet()) {//add inter fragment bonds to cloned fragments
			for (Bond bond : fragToInterFragmentBond.get(originalFragment)) {
				interFragmentBondsToClone.add(bond);
			}
		}
		for (Bond bond : interFragmentBondsToClone) {
			Atom originalFromAtom = bond.getFromAtom();
			Atom originalToAtom = bond.getToAtom();
			Fragment originalFragment1 = originalFromAtom.getFrag();
			Fragment originalFragment2 = originalToAtom.getFrag();
			if (!oldNewFragmentMapping.containsKey(originalFragment1) || (!oldNewFragmentMapping.containsKey(originalFragment2))){
				throw new StructureBuildingException("An element that was clone contained a bond that went outside the scope of the cloning");
			}
			Fragment newFragment1 = oldNewFragmentMapping.get(originalFragment1);
			Fragment newFragment2 = oldNewFragmentMapping.get(originalFragment2);
			Atom fromAtom = newFragment1.getAtomList().get(originalFragment1.getAtomList().indexOf(originalFromAtom));
			Atom toAtom = newFragment2.getAtomList().get(originalFragment2.getAtomList().indexOf(originalToAtom));
			createBond(fromAtom, toAtom, bond.getOrder());
		}
		return clone;
	}

	/**
	 * Takes a terminalAtom on one atom and an atom on another fragment. The properties (charge/element) of the second atom and all of its neighbours are
	 * moved onto the terminalAtom
	 * @param terminalAtom
	 * @param atomThatReplacesTerminal
	 * @throws StructureBuildingException
	 */
	void replaceTerminalAtomWithFragment(Atom terminalAtom, Atom atomThatReplacesTerminal) throws StructureBuildingException {
		Fragment parentFrag = terminalAtom.getFrag();
		Fragment childFrag = atomThatReplacesTerminal.getFrag();
		if (parentFrag == childFrag){
			throw new StructureBuildingException("Replacing atom and terminal should be different fragments");
		}
		List<Atom> atomNeighbours = terminalAtom.getAtomNeighbours();
		if (atomNeighbours.size() > 1){
			throw new StructureBuildingException("Atom to be attached to fragment should only have one bond");
		}
		terminalAtom.setElement(atomThatReplacesTerminal.getElement());
		terminalAtom.setCharge(atomThatReplacesTerminal.getCharge());
		terminalAtom.clearLocants();
		for (String locant : atomThatReplacesTerminal.getLocants()) {
			terminalAtom.addLocant(locant);
		}

		List<Atom> neighbours = atomThatReplacesTerminal.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			createBond(terminalAtom, neighbour, childFrag.findBond(atomThatReplacesTerminal, neighbour).getOrder());
		}
		removeAtomAndAssociatedBonds(atomThatReplacesTerminal);
		incorporateFragment(childFrag, parentFrag);
	}

	/**
	 * Checks if this bond is an inter fragment bond and if it is removes it
	 * @param bond
	 */
	private void removeInterFragmentBondIfPresent(Bond bond) {
		if (bondPile.remove(bond)){
			fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).remove(bond);
			fragToInterFragmentBond.get(bond.getToAtom().getFrag()).remove(bond);
		}
	}
	
	/**
	 * Adds a bond to the inter fragment bond list and fragment to inter-fragment bond mappings
	 * @param bond
	 */
	private void addInterFragmentBond(Bond bond) {
		bondPile.add(bond);
		fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).add(bond);
		fragToInterFragmentBond.get(bond.getToAtom().getFrag()).add(bond);
	}

	/**
	 * Gets a set of the inter fragment bonds a fragment is involved in
	 * @param frag
	 * @return set of inter fragment bonds
	 */
	Set<Bond> getInterFragmentBonds(Fragment frag) {
		return fragToInterFragmentBond.get(frag);
	}

	/**
	 * Create a new Atom of the given element belonging to the given fragment
	 * @param elementSymbol
	 * @param frag
	 * @return Atom
	 * @throws StructureBuildingException
	 */
	Atom createAtom(String elementSymbol, Fragment frag) throws StructureBuildingException {
		Atom a = new Atom(idManager.getNextID(), elementSymbol, frag);
		frag.addAtom(a);
		return a;
	}
	
	/**
	 * Create a new bond between two atoms.
	 * The bond is associated with these atoms.
	 * It is also listed as an inter-fragment bond or associated with a fragment
	 * @param fromAtom
	 * @param toAtom
	 * @param bondOrder
	 * @return Bond
	 */
	Bond createBond(Atom fromAtom, Atom toAtom, int bondOrder) {
		Bond b = new Bond(fromAtom, toAtom, bondOrder);
		fromAtom.addBond(b);
		toAtom.addBond(b);
		if (fromAtom.getFrag() == toAtom.getFrag()){
			fromAtom.getFrag().addBond(b);
		}
		else{
			addInterFragmentBond(b);
		}
		return b;
	}
	
	void removeAtomAndAssociatedBonds(Atom atom){
		ArrayList<Bond> bondsToBeRemoved=new ArrayList<Bond>(atom.getBonds());
		for (Bond bond : bondsToBeRemoved) {
			removeBond(bond);
		}
		atom.getFrag().removeAtom(atom);
	}
	
	void removeBond(Bond bond){
		bond.getFromAtom().getFrag().removeBond(bond);
		bond.getFromAtom().removeBond(bond);
		bond.getToAtom().removeBond(bond);
		removeInterFragmentBondIfPresent(bond);
	}
}