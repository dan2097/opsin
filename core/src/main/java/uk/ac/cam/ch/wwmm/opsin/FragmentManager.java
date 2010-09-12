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
	private final Map<Fragment,LinkedHashSet<Bond>> fragToInterFragmentBond;
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
		fragToInterFragmentBond = new HashMap<Fragment, LinkedHashSet<Bond>>();
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
	 * Also imports outAtoms/inAtoms and functionalAtoms
	 * Reassigns inter fragment bonds of the parent fragment as either intra fragment bonds
	 * of the parent fragment or as inter fragment bonds of the parent fragment
	 *
	 * The original fragment still maintains its original atomList/bondList/interFragmentBondList which is necessary for stereochemistry handling
	 *
	 * @param childFrag The fragment to be incorporated
	 * @param parentFrag The parent fragment
     * @throws StructureBuildingException
	 */
	void incorporateFragment(Fragment childFrag, Fragment parentFrag) throws StructureBuildingException {
		for(Atom atom : childFrag.getAtomList()) {
			parentFrag.addAtom(atom);
		}
		for(Bond bond : childFrag.getBondSet()) {
			parentFrag.addBond(bond);
		}
		parentFrag.addOutAtoms(childFrag.getOutAtoms());
		parentFrag.addInAtoms(childFrag.getInAtoms());
		parentFrag.addFunctionalAtoms(childFrag.getFunctionalAtoms());

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
     * @throws StructureBuildingException
	 */
	void incorporateFragment(Fragment childFrag, int fromID, Fragment parentFrag, int toID, int bondOrder) throws StructureBuildingException {
		incorporateFragment(childFrag, parentFrag);
		createBond(parentFrag.getAtomByIDOrThrow(fromID), parentFrag.getAtomByIDOrThrow(toID), bondOrder);
	}

	/** Converts an atom in a fragment to a different atomic symbol.
	 * Charged atoms can also be specified using a SMILES formula eg. [N+]
	 *
	 * @param a The atom to change to a heteroatom
	 * @param smiles The SMILES for one atom
	 * @param assignLocant Whether a locant should be assigned to the heteroatom if the locant is not used elsewhere
	 * @throws StructureBuildingException if the atom could not be found
	 */
	void makeHeteroatom(Atom a, String smiles, boolean assignLocant) throws StructureBuildingException {
		String elementSymbol;
		if(smiles.startsWith("[")) {
			Atom heteroAtom = sBuilder.build(smiles, this).getFirstAtom();
			elementSymbol =heteroAtom.getElement();
			int charge = heteroAtom.getCharge();
			if (heteroAtom.getCharge()!=0){
				Integer defaultValency = ValencyChecker.getDefaultValency(elementSymbol);
				if (defaultValency !=null){
					//calculate how many protons have been added/removed as compared to the element's standard valency
					int specifiedValency =  heteroAtom.getMinimumValency() !=null ?
							heteroAtom.getMinimumValency() : ValencyChecker.getPossibleValencies(elementSymbol, charge)[0];
					int hydrogensAdded = specifiedValency - defaultValency;
					a.addChargeAndProtons(charge , hydrogensAdded);
				}
				else{
					a.addChargeAndProtons(charge , 0);
				}
			}
		}
		else{
			elementSymbol = smiles;
		}
		a.setElement(elementSymbol);
		a.removeElementSymbolLocants();
		if (assignLocant && a.getFrag().getAtomByLocant(elementSymbol) ==null){//if none of that element currently present add element symbol locant
			a.addLocant(elementSymbol);
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
     * @throws StructureBuildingException
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
			FragmentTools.convertSpareValenciesToDoubleBonds(f);
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
		fragToInterFragmentBond.put(frag, new LinkedHashSet<Bond>());
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
	Fragment copyFragment(Fragment originalFragment) throws StructureBuildingException {
		return copyAndRelabelFragment(originalFragment, null);
	}


	/**
	 * Creates a copy of a fragment by copying data
	 * labels the atoms using new ids from the idManager
	 * @param originalFragment
	 * @param stringToAddToAllLocants: typically used to append primes to all locants, can be null
	 * @return the clone of the fragment
	 * @throws StructureBuildingException
	 */
	Fragment copyAndRelabelFragment(Fragment originalFragment, String stringToAddToAllLocants) throws StructureBuildingException {
		Fragment newFragment =new Fragment(originalFragment.getType(), originalFragment.getSubType());
		HashMap<Atom, Atom> oldToNewAtomMap = new HashMap<Atom, Atom>();//maps old Atom to new Atom
		List<Atom> atomList =originalFragment.getAtomList();
		List<OutAtom> outAtoms =originalFragment.getOutAtoms();
		List<FunctionalAtom> functionalAtoms =originalFragment.getFunctionalAtoms();
		List<InAtom> inAtoms =originalFragment.getInAtoms();
		Atom defaultInAtom =originalFragment.getDefaultInAtom();
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
			newAtom.setLambdaConventionValency(atom.getLambdaConventionValency());
			//outValency is derived from the outAtoms so is automatically cloned
			newAtom.setAtomIsInACycle(atom.getAtomIsInACycle());
			newAtom.setType(atom.getType());//may be different from fragment type if the original atom was formerly in a suffix
			newAtom.setMinimumValency(atom.getMinimumValency());
			newFragment.addAtom(newAtom);
			oldToNewAtomMap.put(atom, newAtom);
		}
		for (Atom atom : atomList) {
			if (atom.getAtomParity() != null){
				Atom[] oldAtomRefs4 = atom.getAtomParity().getAtomRefs4();
				Atom[] newAtomRefs4 = new Atom[4];
				for (int i = 0; i < oldAtomRefs4.length; i++) {
					Atom oldAtom = oldAtomRefs4[i];
					if (oldAtom.equals(AtomParity.hydrogen)){
						newAtomRefs4[i] = AtomParity.hydrogen;
					}
					else{
						newAtomRefs4[i] = oldToNewAtomMap.get(oldAtom);
					}
				}
				AtomParity newAtomParity =new AtomParity(newAtomRefs4, atom.getAtomParity().getParity());
				oldToNewAtomMap.get(atom).setAtomParity(newAtomParity);
			}
			if (atom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
				Set<Atom> oldAtoms = atom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
				Set<Atom> newAtoms = new HashSet<Atom>();
				for (Atom oldAtom : oldAtoms) {
					newAtoms.add(oldToNewAtomMap.get(oldAtom));
				}
				oldToNewAtomMap.get(atom).setProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT, newAtoms);
			}
			if (atom.getProperty(Atom.SMILES_HYDROGEN_COUNT)!=null){
				oldToNewAtomMap.get(atom).setProperty(Atom.SMILES_HYDROGEN_COUNT, atom.getProperty(Atom.SMILES_HYDROGEN_COUNT));
			}
		}
        for (OutAtom outAtom : outAtoms) {
            newFragment.addOutAtom(oldToNewAtomMap.get(outAtom.getAtom()), outAtom.getValency(), outAtom.isSetExplicitly());
        }
		for (FunctionalAtom functionalAtom : functionalAtoms) {
			newFragment.addFunctionalAtom(oldToNewAtomMap.get(functionalAtom.getAtom()));
		}
        for (InAtom inAtom : inAtoms) {
            newFragment.addInAtom(oldToNewAtomMap.get(inAtom.getAtom()), inAtom.getValency());
        }
		newFragment.setDefaultInAtom(oldToNewAtomMap.get(defaultInAtom));
		Set<Bond> bondSet =originalFragment.getBondSet();
		for (Bond bond : bondSet) {
			Bond newBond = createBond(oldToNewAtomMap.get(bond.getFromAtom()), oldToNewAtomMap.get(bond.getToAtom()), bond.getOrder());
			newBond.setSmilesStereochemistry(bond.getSmilesStereochemistry());
			if (bond.getBondStereoElement() != null){
				newBond.setBondStereoElement(new Element(bond.getBondStereoElement()));
			}
		}
		List<Atom> indicatedHydrogenAtoms = originalFragment.getIndicatedHydrogen();
		for (Atom atom : indicatedHydrogenAtoms) {
			newFragment.addIndicatedHydrogen(oldToNewAtomMap.get(atom));
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
		List<Element> originalGroups = XOMTools.getDescendantElementsWithTagName(elementToBeCloned, XmlDeclarations.GROUP_EL);
		List<Element> clonedGroups = XOMTools.getDescendantElementsWithTagName(clone,  XmlDeclarations.GROUP_EL);
		HashMap<Fragment,Fragment> oldNewFragmentMapping  =new HashMap<Fragment, Fragment>();
		for (int i = 0; i < originalGroups.size(); i++) {
			Fragment originalFragment =state.xmlFragmentMap.get(originalGroups.get(i));
			Fragment newFragment = copyAndRelabelFragment(originalFragment, stringToAddToAllLocants);
			oldNewFragmentMapping.put(originalFragment, newFragment);
			state.xmlFragmentMap.put(clonedGroups.get(i), newFragment);
			List<Fragment> originalSuffixes =state.xmlSuffixMap.get(originalGroups.get(i));
			List<Fragment> newSuffixFragments =new ArrayList<Fragment>();
			for (Fragment suffix : originalSuffixes) {
				newSuffixFragments.add(state.fragManager.copyFragment(suffix));
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
				throw new StructureBuildingException("An element that was a clone contained a bond that went outside the scope of the cloning");
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
		Fragment parentFrag = terminalAtom.getFrag();//TODO  use replaceAtomWithAnotherAtomPreservingConnectivity instead of this function????
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
	 * Takes an atom, removes it and bonds everything that was bonded to it to the replacementAtom with the original bond orders.
	 * Non element symbol locants are copied to the replacement atom
	 * @param atomToBeReplaced
	 * @param replacementAtom
	 */
	void replaceAtomWithAnotherAtomPreservingConnectivity(Atom atomToBeReplaced, Atom replacementAtom) {
		atomToBeReplaced.removeElementSymbolLocants();
		List<String> locants = atomToBeReplaced.getLocants();
		for (int i = locants.size() -1; i >=0; i--) {
			String locant = locants.get(i);
			atomToBeReplaced.removeLocant(locant);
			replacementAtom.addLocant(locant);
		}
		Fragment replacedAtomsFragment = atomToBeReplaced.getFrag();
		List<Atom> neighbours = atomToBeReplaced.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			createBond(replacementAtom, neighbour, replacedAtomsFragment.findBond(atomToBeReplaced, neighbour).getOrder());
		}
		removeAtomAndAssociatedBonds(atomToBeReplaced);
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