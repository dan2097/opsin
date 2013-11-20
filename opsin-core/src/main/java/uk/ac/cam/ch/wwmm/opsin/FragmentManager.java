package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import nu.xom.Element;


/** Holds the Fragments during the construction of the molecule,
 * handles the building of new fragments and handles the creation/deletion of atoms/bonds
 *
 * @author ptc24
 * @author dl387
 *
 */
class FragmentManager {

	/** A mapping between fragments and inter fragment bonds */
	private final Map<Fragment,LinkedHashSet<Bond>> fragToInterFragmentBond = new LinkedHashMap<Fragment, LinkedHashSet<Bond>>();
	
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
		Fragment newFrag = sBuilder.build(smiles, type, subType, labelMapping);
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
	 * @throws StructureBuildingException 
	 */
	Fragment getUnifiedFragment() throws StructureBuildingException {
		Fragment uniFrag = new Fragment();
		for (Entry<Fragment, LinkedHashSet<Bond>> entry : fragToInterFragmentBond.entrySet()) {
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
		String elementSymbol =heteroAtom.getElement();
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
		a.setElement(elementSymbol);
		a.removeElementSymbolLocants();
		if (assignLocant){
			String primes ="";
			while (a.getFrag().getAtomByLocant(elementSymbol+primes)!=null){//if element symbol already assigned, add a prime and try again
				primes+="'";
			}
			a.addLocant(elementSymbol +primes);
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
	 * @throws StructureBuildingException
	 */
	Fragment copyAndRelabelFragment(Fragment originalFragment, int primesToAdd) throws StructureBuildingException {
		Fragment newFragment =new Fragment(originalFragment.getType(), originalFragment.getSubType());
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
				Set<Atom> newAtoms = new HashSet<Atom>();
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
		Element clone = new Element(elementToBeCloned);
		List<Element> originalGroups = XOMTools.getDescendantElementsWithTagName(elementToBeCloned, XmlDeclarations.GROUP_EL);
		List<Element> clonedGroups = XOMTools.getDescendantElementsWithTagName(clone,  XmlDeclarations.GROUP_EL);
		HashMap<Fragment,Fragment> oldNewFragmentMapping  =new LinkedHashMap<Fragment, Fragment>();
		for (int i = 0; i < originalGroups.size(); i++) {
			Fragment originalFragment =state.xmlFragmentMap.get(originalGroups.get(i));
			Fragment newFragment = copyAndRelabelFragment(originalFragment, primesToAdd);
			oldNewFragmentMapping.put(originalFragment, newFragment);
			state.xmlFragmentMap.put(clonedGroups.get(i), newFragment);
			List<Fragment> originalSuffixes =state.xmlSuffixMap.get(originalGroups.get(i));
			List<Fragment> newSuffixFragments =new ArrayList<Fragment>();
			for (Fragment suffix : originalSuffixes) {
				newSuffixFragments.add(state.fragManager.copyFragment(suffix));
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
		List<Bond> bondsToBeRemoved = new ArrayList<Bond>(atom.getBonds());
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