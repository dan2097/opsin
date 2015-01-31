package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/** Holds the Fragments during the construction of the molecule,
 * handles the building of new fragments and handles the creation/deletion of atoms/bonds
 *
 * @author ptc24
 * @author dl387
 *
 */
class FragmentManager {

	/** A mapping between fragments and inter fragment bonds */
	private final Map<Fragment,Set<Bond>> fragToInterFragmentBond = new LinkedHashMap<Fragment, Set<Bond>>();
	
	/** All of the atom-containing fragments in the molecule */
	private final Set<Fragment> fragments = fragToInterFragmentBond.keySet();
	
	/** A builder for fragments specified as SMILES */
	private final SMILESFragmentBuilder sBuilder;
	
	/** A source of unique integers */
	private final IDManager idManager;

	/** Sets up a new Fragment manager, containing no fragments.
	 *
	 * @param sBuilder A SMILESFragmentBuilder - dependency injection.
	 * @param idManager An IDManager.
	 */
	FragmentManager(SMILESFragmentBuilder sBuilder, IDManager idManager) {
		if (sBuilder == null || idManager == null ){
			throw new IllegalArgumentException("FragmentManager was parsed a null object in its constructor!");
		}
		this.sBuilder = sBuilder;
		this.idManager = idManager;
	}

	/** Builds a fragment, based on an SMILES string
	 * The fragment will not correspond to a token
	 *
	 * @param smiles The fragment to build
	 * @return The built fragment
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles) throws StructureBuildingException {
		return buildSMILES(smiles, "", "");
	}
	
	/** Builds a fragment, based on an SMILES string
	 * The fragment will not correspond to a token
	 * 
	 * @param smiles
	 * @param type
	 * @param labelMapping
	 * @return
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles, String type, String labelMapping) throws StructureBuildingException {
		Fragment newFrag = sBuilder.build(smiles, type, labelMapping);
		addFragment(newFrag);
		return newFrag;
	}

	/** Builds a fragment, based on an SMILES string
	 * The fragment will correspond to the given tokenEl
	 * 
	 * @param smiles The fragment to build
	 * @param tokenEl The corresponding tokenEl
	 * @param labelMapping How to label the fragment
	 * @return The built fragment
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles, Element tokenEl, String labelMapping) throws StructureBuildingException {
		Fragment newFrag = sBuilder.build(smiles, tokenEl, labelMapping);
		addFragment(newFrag);
		return newFrag;
	}

	/**Creates a new fragment, containing all of the atoms and bonds
	 * of all of the other fragments - i.e. the whole molecule. This updates
	 * which fragments the atoms think they are in to the new super fragment
	 * but does not change the original fragments.
	 * Hence the original fragments remain associated with their atoms
	 * Atoms and Bonds are not copied.
	 *
	 * @return The unified fragment
	 */
	Fragment getUnifiedFragment() {
		Fragment uniFrag = new Fragment("");
		for (Entry<Fragment, Set<Bond>> entry : fragToInterFragmentBond.entrySet()) {
			Fragment f = entry.getKey();
			Set<Bond> interFragmentBonds = entry.getValue();
			for(Atom atom : f.getAtomList()) {
				uniFrag.addAtom(atom);
			}
			for(Bond bond : f.getBondSet()) {
				uniFrag.addBond(bond);
			}
			uniFrag.incorporateOutAtoms(f);
			uniFrag.incorporateFunctionalAtoms(f);

			for (Bond interFragmentBond : interFragmentBonds) {
				uniFrag.addBond(interFragmentBond);
			}
		}
		addFragment(uniFrag);
		return uniFrag;
	}

	/** Incorporates a fragment, usually a suffix, into a parent fragment
	 * This does:
	 * Imports all of the atoms and bonds from another fragment into this one.
	 * Also imports outAtoms and functionalAtoms
	 * Reassigns inter-fragment bonds of the child fragment as either intra-fragment bonds
	 * of the parent fragment or as inter-fragment bonds of the parent fragment
	 *
	 * The original fragment still maintains its original atomList/bondList
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
		parentFrag.incorporateOutAtoms(childFrag);
		parentFrag.incorporateFunctionalAtoms(childFrag);

		Set<Bond> interFragmentBonds = fragToInterFragmentBond.get(childFrag);
		if (interFragmentBonds == null){
			throw new StructureBuildingException("Fragment not registered with this FragmentManager!");
		}
		for (Bond bond : interFragmentBonds) {//reassign inter-fragment bonds of child
			if (bond.getFromAtom().getFrag() == parentFrag && bond.getToAtom().getFrag() == parentFrag){
				//bond is now enclosed within parentFrag so make it an intra-fragment bond
				//and remove it from the inter-fragment set of the parentFrag
				parentFrag.addBond(bond);
				fragToInterFragmentBond.get(parentFrag).remove(bond);
			}
			else{
				//bond was an inter-fragment bond between the childFrag and another frag
				//It is now between the parentFrag and another frag
				addInterFragmentBond(bond);
			}
		}
		fragToInterFragmentBond.remove(childFrag);
	}
	
	/** Incorporates a fragment, usually a suffix, into a parent fragment, creating a bond between them.
	 *
	 * @param childFrag The fragment to be incorporated
	 * @param fromAtom An atom on that fragment
	 * @param parentFrag The parent fragment
	 * @param toAtom An atom on that fragment
	 * @param bondOrder The order of the joining bond
     * @throws StructureBuildingException
	 */
	void incorporateFragment(Fragment childFrag, Atom fromAtom, Fragment parentFrag, Atom toAtom, int bondOrder) throws StructureBuildingException {
		if (!fromAtom.getFrag().equals(childFrag)){
			throw new StructureBuildingException("OPSIN Bug: fromAtom was not associated with childFrag!");
		}
		if (!toAtom.getFrag().equals(parentFrag)){
			throw new StructureBuildingException("OPSIN Bug: toAtom was not associated with parentFrag!");
		}
		incorporateFragment(childFrag, parentFrag);
		createBond(fromAtom, toAtom, bondOrder);
	}

	/** Converts an atom in a fragment to a different atomic symbol described by a SMILES string
	 * Charged atoms can also be specified eg. [NH4+]
	 *
	 * @param a The atom to change to a heteroatom
	 * @param smiles The SMILES for one atom
	 * @throws StructureBuildingException
	 */
	void replaceAtomWithSmiles(Atom a, String smiles) throws StructureBuildingException {
		replaceAtomWithAtom(a, getHeteroatom(smiles), false);
	}

	/**
	 * Converts the smiles for a heteroatom to an atom
	 * @param smiles
	 * @return
	 * @throws StructureBuildingException
	 */
	Atom getHeteroatom(String smiles) throws StructureBuildingException {
		Fragment heteroAtomFrag = sBuilder.build(smiles);
		if (heteroAtomFrag.getAtomCount() != 1){
			throw new StructureBuildingException("Heteroatom smiles described a fragment with multiple SMILES!");
		}
		return heteroAtomFrag.getFirstAtom();
	}
	
	/** Uses the information given in the given heteroatom to change the atomic symbol
	 * and charge of the given atom
	 *
	 * @param a The atom to change to a heteroatom
	 * @param heteroAtom The atom to copy element/charge properties from
	 * @param assignLocant Whether a locant should be assigned to the heteroatom if the locant is not used elsewhere
	 * @throws StructureBuildingException if a charge disagreement occurs
	 */
	void replaceAtomWithAtom(Atom a, Atom heteroAtom, boolean assignLocant) throws StructureBuildingException {
		ChemEl chemEl =heteroAtom.getElement();
		int replacementCharge =heteroAtom.getCharge();
		if (replacementCharge!=0){
			if (a.getCharge()==0){
				a.addChargeAndProtons(replacementCharge, heteroAtom.getProtonsExplicitlyAddedOrRemoved());
			}
			else if (a.getCharge()==replacementCharge){
				a.setProtonsExplicitlyAddedOrRemoved(heteroAtom.getProtonsExplicitlyAddedOrRemoved());
			}
			else{
				throw new StructureBuildingException("Charge conflict between replacement term and atom to be replaced");
			}
		}
		a.setElement(chemEl);
		a.removeElementSymbolLocants();
		if (assignLocant){
			String primes = "";
			while (a.getFrag().getAtomByLocant(chemEl.toString() + primes) != null){//if element symbol already assigned, add a prime and try again
				primes += "'";
			}
			a.addLocant(chemEl.toString() + primes);
		}
	}

	/** Gets an atom, given an id number
	 * Use this if you don't know what fragment the atom is in
	 * @param id The id of the atom
	 * @return The atom, or null if no such atom exists.
	 */
	Atom getAtomByID(int id) {
		for(Fragment f : fragments) {
			Atom a = f.getAtomByID(id);
			if(a != null) {
				return a;
			}
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
		if(a == null) {
			throw new StructureBuildingException("Couldn't get atom by id");
		}
		return a;
	}

	/**Turns all of the spare valencies in the fragments into double bonds.
	 *
	 * @throws StructureBuildingException
	 */
	void convertSpareValenciesToDoubleBonds() throws StructureBuildingException {
		for(Fragment f : fragments) {
			FragmentTools.convertSpareValenciesToDoubleBonds(f);
		}
	}

	/**
	 * Checks valencies are all chemically reasonable. An exception is thrown if any are not
	 * @throws StructureBuildingException
	 */
	void checkValencies() throws StructureBuildingException {
		for(Fragment f : fragments) {
			f.checkValencies();
		}
	}

	Set<Fragment> getFragments() {
		return Collections.unmodifiableSet(fragments);
	}

	/**
	 * Registers a fragment
	 * @param frag
	 */
	private void addFragment(Fragment frag)  {
		fragToInterFragmentBond.put(frag, new LinkedHashSet<Bond>());
	}

	/**
	 * Removes a fragment
	 * Any inter-fragment bonds of this fragment are removed from the fragments it was connected to
	 * Throws an exception if fragment wasn't present
	 * @param frag
	 * @throws StructureBuildingException
	 */
	void removeFragment(Fragment frag) throws StructureBuildingException {
		Set<Bond> interFragmentBondsInvolvingFragmentSet = fragToInterFragmentBond.get(frag);
		if (interFragmentBondsInvolvingFragmentSet == null) {
			throw new StructureBuildingException("Fragment not registered with this FragmentManager!");
		}
		List<Bond> interFragmentBondsInvolvingFragment = new ArrayList<Bond>(interFragmentBondsInvolvingFragmentSet);
		for (Bond bond : interFragmentBondsInvolvingFragment) {
			if (bond.getFromAtom().getFrag() == frag){
				fragToInterFragmentBond.get(bond.getToAtom().getFrag()).remove(bond);
			}
			else{
				fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).remove(bond);
			}
		}
		fragToInterFragmentBond.remove(frag);
	}

	int getOverallCharge() {
		int totalCharge = 0;
		for (Fragment frag : fragments) {
			totalCharge += frag.getCharge();
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
		return copyAndRelabelFragment(originalFragment, 0);
	}

	/**
	 * Creates a copy of a fragment by copying data
	 * labels the atoms using new ids from the idManager
	 * @param originalFragment
	 * @param primesToAdd: The minimum number of primes to add to the cloned atoms. More primes will be added if necessary to keep the locants unique e.g. N in the presence of N' becomes N'' when this is 1
	 * @return the clone of the fragment
	 */
	Fragment copyAndRelabelFragment(Fragment originalFragment, int primesToAdd) {
		Element tokenEl = new TokenEl("");
		tokenEl.addAttribute(TYPE_ATR, originalFragment.getType());
		tokenEl.addAttribute(SUBTYPE_ATR, originalFragment.getSubType());
		Fragment newFragment = new Fragment(tokenEl);
		HashMap<Atom, Atom> oldToNewAtomMap = new HashMap<Atom, Atom>();//maps old Atom to new Atom
		List<Atom> atomList =originalFragment.getAtomList();
		Atom defaultInAtom =originalFragment.getDefaultInAtom();
		for (Atom atom : atomList) {
			int id = idManager.getNextID();
			ArrayList<String> newLocants = new ArrayList<String>(atom.getLocants());
			if (primesToAdd !=0){
				for (int i = 0; i < newLocants.size(); i++) {
					String currentLocant = newLocants.get(i);
					int currentPrimes = StringTools.countTerminalPrimes(currentLocant);
					String locantSansPrimes = currentLocant.substring(0, currentLocant.length()-currentPrimes);
					int highestNumberOfPrimesWithThisLocant = currentPrimes;
					while (originalFragment.getAtomByLocant(locantSansPrimes + StringTools.multiplyString("'", highestNumberOfPrimesWithThisLocant +1 ))!=null){
						highestNumberOfPrimesWithThisLocant++;
					}
					newLocants.set(i, locantSansPrimes + StringTools.multiplyString("'", ((highestNumberOfPrimesWithThisLocant +1)*primesToAdd) + currentPrimes));
				}
			}
			Atom newAtom =new Atom(id, atom.getElement(), newFragment);
			for (String newLocant : newLocants) {
				newAtom.addLocant(newLocant);
			}
			newAtom.setCharge(atom.getCharge());
			newAtom.setIsotope(atom.getIsotope());
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
					else if (oldAtom.equals(AtomParity.deoxyHydrogen)){
						newAtomRefs4[i] = AtomParity.deoxyHydrogen;
					}
					else{
						newAtomRefs4[i] = oldToNewAtomMap.get(oldAtom);
					}
				}
				AtomParity newAtomParity =new AtomParity(newAtomRefs4, atom.getAtomParity().getParity());
				oldToNewAtomMap.get(atom).setAtomParity(newAtomParity);
			}
			Set<Atom> oldAmbiguousElementAssignmentAtoms = atom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
			if (oldAmbiguousElementAssignmentAtoms!=null){
				Set<Atom> newAtoms = new LinkedHashSet<Atom>();
				for (Atom oldAtom : oldAmbiguousElementAssignmentAtoms) {
					newAtoms.add(oldToNewAtomMap.get(oldAtom));
				}
				oldToNewAtomMap.get(atom).setProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT, newAtoms);
			}
			Integer smilesHydrogenCount = atom.getProperty(Atom.SMILES_HYDROGEN_COUNT);
			if (smilesHydrogenCount!=null){
				oldToNewAtomMap.get(atom).setProperty(Atom.SMILES_HYDROGEN_COUNT, smilesHydrogenCount);
			}
			Integer oxidationNumber = atom.getProperty(Atom.OXIDATION_NUMBER);
			if (oxidationNumber!=null){
				oldToNewAtomMap.get(atom).setProperty(Atom.OXIDATION_NUMBER, oxidationNumber);
			}
			Boolean isAldehyde = atom.getProperty(Atom.ISALDEHYDE);
			if (isAldehyde!=null){
				oldToNewAtomMap.get(atom).setProperty(Atom.ISALDEHYDE, isAldehyde);
			}
			Boolean isAnomeric = atom.getProperty(Atom.ISANOMERIC);
			if (isAnomeric!=null){
				oldToNewAtomMap.get(atom).setProperty(Atom.ISANOMERIC, isAnomeric);
			}
		}
		for (int i = 0, l = originalFragment.getOutAtomCount(); i < l; i++) {
			OutAtom outAtom = originalFragment.getOutAtom(i);
            newFragment.addOutAtom(oldToNewAtomMap.get(outAtom.getAtom()), outAtom.getValency(), outAtom.isSetExplicitly());
            if (outAtom.getLocant() !=null){
            	newFragment.getOutAtom(newFragment.getOutAtomCount() -1).setLocant(outAtom.getLocant() + StringTools.multiplyString("'", primesToAdd) );
            }
        }
		for (int i = 0, l = originalFragment.getFunctionalAtomCount(); i < l; i++) {
			FunctionalAtom functionalAtom = originalFragment.getFunctionalAtom(i);
			newFragment.addFunctionalAtom(oldToNewAtomMap.get(functionalAtom.getAtom()));
		}
		newFragment.setDefaultInAtom(oldToNewAtomMap.get(defaultInAtom));
		Set<Bond> bondSet =originalFragment.getBondSet();
		for (Bond bond : bondSet) {
			Bond newBond = createBond(oldToNewAtomMap.get(bond.getFromAtom()), oldToNewAtomMap.get(bond.getToAtom()), bond.getOrder());
			newBond.setSmilesStereochemistry(bond.getSmilesStereochemistry());
			if (bond.getBondStereo() != null){
				Atom[] oldAtomRefs4 = bond.getBondStereo().getAtomRefs4();
				Atom[] newAtomRefs4 = new Atom[4];
				for (int i = 0; i < oldAtomRefs4.length; i++) {
					newAtomRefs4[i] = oldToNewAtomMap.get(oldAtomRefs4[i]);
				}
				newBond.setBondStereoElement(newAtomRefs4, bond.getBondStereo().getBondStereoValue());
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
		return cloneElement(state, elementToBeCloned, 0);
	}

	/**
	 * Takes an element and produces a copy of it. Groups and suffixes are copied so that the new element
	 * has its own group and suffix fragments
	 * @param elementToBeCloned
	 * @param state The current buildstate
	 * @param primesToAdd: The minimum number of primes to add to the cloned atoms. More primes will be added if necessary to keep the locants unique e.g. N in the presence of N' becomes N'' when this is 1
	 * @return
	 * @throws StructureBuildingException
	 */
	Element cloneElement(BuildState state, Element elementToBeCloned, int primesToAdd) throws StructureBuildingException {
		Element clone = elementToBeCloned.copy();
		List<Element> originalGroups = OpsinTools.getDescendantElementsWithTagName(elementToBeCloned, XmlDeclarations.GROUP_EL);
		List<Element> clonedGroups = OpsinTools.getDescendantElementsWithTagName(clone,  XmlDeclarations.GROUP_EL);
		HashMap<Fragment,Fragment> oldNewFragmentMapping  =new LinkedHashMap<Fragment, Fragment>();
		for (int i = 0; i < originalGroups.size(); i++) {
			Fragment originalFragment = originalGroups.get(i).getFrag();
			Fragment newFragment = copyAndRelabelFragment(originalFragment, primesToAdd);
			oldNewFragmentMapping.put(originalFragment, newFragment);
			newFragment.setTokenEl(clonedGroups.get(i));
			clonedGroups.get(i).setFrag(newFragment);
			List<Fragment> originalSuffixes =state.xmlSuffixMap.get(originalGroups.get(i));
			List<Fragment> newSuffixFragments =new ArrayList<Fragment>();
			for (Fragment suffix : originalSuffixes) {
				newSuffixFragments.add(copyFragment(suffix));
			}
			state.xmlSuffixMap.put(clonedGroups.get(i), newSuffixFragments);
		}
		Set<Bond> interFragmentBondsToClone = new LinkedHashSet<Bond>();
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
	 * Takes an atom, removes it and bonds everything that was bonded to it to the replacementAtom with the original bond orders.
	 * Non element symbol locants are copied to the replacement atom
	 * @param atomToBeReplaced
	 * @param replacementAtom
	 */
	void replaceAtomWithAnotherAtomPreservingConnectivity(Atom atomToBeReplaced, Atom replacementAtom) {
		atomToBeReplaced.removeElementSymbolLocants();
		List<String> locants = new ArrayList<String>(atomToBeReplaced.getLocants());
		for (String locant : locants) {
			atomToBeReplaced.removeLocant(locant);
			replacementAtom.addLocant(locant);
		}
		List<Bond> bonds = atomToBeReplaced.getBonds();
		for (Bond bond : bonds) {
			Atom connectedAtom = bond.getOtherAtom(atomToBeReplaced);
			if (connectedAtom.getAtomParity() != null){
				Atom[] atomRefs4 = connectedAtom.getAtomParity().getAtomRefs4();
				for (int i = 0 ; i < 4; i++) {
					if (atomRefs4[i] == atomToBeReplaced){
						atomRefs4[i] = replacementAtom;
						break;
					}
				}
			}
			if (bond.getBondStereo() != null){
				Atom[] atomRefs4 = bond.getBondStereo().getAtomRefs4();
				for (int i = 0 ; i < 4; i++) {
					if (atomRefs4[i] == atomToBeReplaced){
						atomRefs4[i] = replacementAtom;
						break;
					}
				}
			}
			createBond(replacementAtom, bond.getOtherAtom(atomToBeReplaced), bond.getOrder());
		}
		removeAtomAndAssociatedBonds(atomToBeReplaced);
	}

	/**
	 * Removes a bond from the inter-fragment bond mappings if it was present
	 * @param bond
	 */
	private void removeInterFragmentBondIfPresent(Bond bond) {
		fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).remove(bond);
		fragToInterFragmentBond.get(bond.getToAtom().getFrag()).remove(bond);
	}
	
	/**
	 * Adds a bond to the fragment to inter-fragment bond mappings
	 * @param bond
	 */
	private void addInterFragmentBond(Bond bond) {
		fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).add(bond);
		fragToInterFragmentBond.get(bond.getToAtom().getFrag()).add(bond);
	}

	/**
	 * Gets an unmodifiable view of the set of the inter-fragment bonds a fragment is involved in
	 * @param frag
	 * @return set of inter fragment bonds
	 */
	Set<Bond> getInterFragmentBonds(Fragment frag) {
		Set<Bond> interFragmentBonds = fragToInterFragmentBond.get(frag);
		if (interFragmentBonds == null) {
			throw new IllegalArgumentException("Fragment not registered with this FragmentManager!");
		}
		return Collections.unmodifiableSet(interFragmentBonds);
	}

	/**
	 * Create a new Atom of the given element belonging to the given fragment
	 * @param chemEl
	 * @param frag
	 * @return Atom
	 */
	Atom createAtom(ChemEl chemEl, Fragment frag) {
		Atom a = new Atom(idManager.getNextID(), chemEl, frag);
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
		List<Bond> bondsToBeRemoved = new ArrayList<Bond>(atom.getBonds());
		for (Bond bond : bondsToBeRemoved) {
			removeBond(bond);
		}
		atom.getFrag().removeAtom(atom);
		Set<Atom> ambiguousElementAssignment = atom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
		if (ambiguousElementAssignment != null){
			ambiguousElementAssignment.remove(atom);
			if (ambiguousElementAssignment.size() == 1){
				ambiguousElementAssignment.iterator().next().setProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT, null);
			}
		}
	}
	
	void removeBond(Bond bond){
		bond.getFromAtom().getFrag().removeBond(bond);
		bond.getFromAtom().removeBond(bond);
		bond.getToAtom().removeBond(bond);
		removeInterFragmentBondIfPresent(bond);
	}
	
	/**
	 * Valency is used to determine the expected number of hydrogen
	 * Hydrogens are then added to bring the number of connections up to the minimum required to satisfy the atom's valency
	 * This allows the valency of the atom to be encoded e.g. phopshane-3 hydrogen, phosphorane-5 hydrogen.
	 * It is also necessary when considering stereochemistry as a hydrogen beats nothing in the CIP rules
	 * @throws StructureBuildingException
	 */
	void makeHydrogensExplicit() throws StructureBuildingException {
		for (Fragment fragment : fragments) {
			if (fragment.getSubType().equals(ELEMENTARYATOM_SUBTYPE_VAL)){//these do not have implicit hydrogen e.g. phosphorus is literally just a phosphorus atom
				continue;
			}
			List<Atom> atomList = fragment.getAtomList();
			for (Atom parentAtom : atomList) {
				int explicitHydrogensToAdd = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(parentAtom);
				for (int i = 0; i < explicitHydrogensToAdd; i++) {
					Atom hydrogen = createAtom(ChemEl.H, fragment);
					createBond(parentAtom, hydrogen, 1);
				}
				if (parentAtom.getAtomParity() != null){
					if (explicitHydrogensToAdd > 1) {
						//Cannot have tetrahedral chirality and more than 2 hydrogens
						parentAtom.setAtomParity(null);//probably caused by deoxy
					}
					else {
						modifyAtomParityToTakeIntoAccountExplicitHydrogen(parentAtom);
					}
				}
			}
		}
	}
	
	private void modifyAtomParityToTakeIntoAccountExplicitHydrogen(Atom atom) throws StructureBuildingException {
		AtomParity atomParity = atom.getAtomParity();
		if (!StereoAnalyser.isPossiblyStereogenic(atom)){
			//no longer a stereoCentre e.g. due to unsaturation
			atom.setAtomParity(null);
		}
		else{
			Atom[] atomRefs4 = atomParity.getAtomRefs4();
			Integer positionOfImplicitHydrogen = null;
			Integer positionOfDeoxyHydrogen = null;
			for (int i = 0; i < atomRefs4.length; i++) {
				Atom a = atomRefs4[i];
				if (a.equals(AtomParity.hydrogen)){
					positionOfImplicitHydrogen = i;
				}
				else if (a.equals(AtomParity.deoxyHydrogen)){
					positionOfDeoxyHydrogen = i;
				}
			}
			if (positionOfImplicitHydrogen != null || positionOfDeoxyHydrogen != null) {
				//atom parity was set in SMILES, the dummy hydrogen atom has now been substituted
				List<Atom> neighbours = atom.getAtomNeighbours();
				for (Atom atomRef : atomRefs4) {
					neighbours.remove(atomRef);
				}
				if (neighbours.size() == 0) {
					throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substituted a hydrogen at stereocentre");
				}
				else if (neighbours.size() == 1 && positionOfDeoxyHydrogen != null) {
					atomRefs4[positionOfDeoxyHydrogen] = neighbours.get(0);
					if (positionOfImplicitHydrogen != null){
						throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substituted a hydrogen at stereocentre");
					}
				}
				else if (neighbours.size() == 1 && positionOfImplicitHydrogen != null) {
					atomRefs4[positionOfImplicitHydrogen] = neighbours.get(0);
				}
				else if (neighbours.size() == 2 && positionOfDeoxyHydrogen != null && positionOfImplicitHydrogen != null) {
					try{
						List<Atom> cipOrderedAtoms = new CipSequenceRules(atom).getNeighbouringAtomsInCIPOrder();
						//higher priority group replaces the former hydroxy groups (deoxyHydrogen)
						if (cipOrderedAtoms.indexOf(neighbours.get(0)) > cipOrderedAtoms.indexOf(neighbours.get(1))) {
							atomRefs4[positionOfDeoxyHydrogen] = neighbours.get(0);
							atomRefs4[positionOfImplicitHydrogen] = neighbours.get(1);
						}
						else{
							atomRefs4[positionOfDeoxyHydrogen] = neighbours.get(1);
							atomRefs4[positionOfImplicitHydrogen] = neighbours.get(0);
						}
					}
					catch (CipOrderingException e){
						//assume ligands equivalent so it makes no difference which is which
						atomRefs4[positionOfDeoxyHydrogen] = neighbours.get(0);
						atomRefs4[positionOfImplicitHydrogen] = neighbours.get(1);
					}
				}
				else{
					throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substituted a hydrogen at stereocentre");
				}
			}
		}
	}
}