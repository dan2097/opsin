package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nu.xom.Element;


/** Holds the Fragments during the construction of the molecule,
 * and handles the building of new fragments.
 *
 * @author ptc24
 *
 */
class FragmentManager {

	/**
	 * Sorts a list of atoms such that their order agrees with the order symbolic locants are typically assigned
	 * @author dl387
	 *
	 */
	static class SortAtomsForElementSymbols implements Comparator<Atom> {

	    public int compare(Atom a, Atom b){
	    	int compare =a.getElement().compareTo(b.getElement());
	    	if (compare !=0){//only bother comparing properly if elements are the same
	    		return compare;
	    	}

	    	if (a.getBonds().size() >b.getBonds().size()){//less bonds is preferred
	    		return 1;
	    	}
	    	if (a.getBonds().size() <b.getBonds().size()){
	    		return -1;
	    	}

	    	int aTotalBondOrder =0;
	    	for (Bond bonda : a.getBonds()) {
				aTotalBondOrder +=bonda.getOrder();
			}

	    	int bTotalBondOrder =0;
	    	for (Bond bondb : b.getBonds()) {
				bTotalBondOrder +=bondb.getOrder();
			}

	    	aTotalBondOrder +=a.getOutValency();//take into account the bond/s which hasn't been added yet which will go to a main fragment
	    	bTotalBondOrder +=b.getOutValency();//take into account the bond/s which hasn't been added yet which will go to a main fragment

	    	if (aTotalBondOrder >bTotalBondOrder){//lower order bonds are preferred
	    		return 1;
	    	}
	    	if (aTotalBondOrder <bTotalBondOrder){
	    		return -1;
	    	}

	    	return 0;
	    }
	}

	/** All of the atom-containing fragments in the molecule */
	private Set<Fragment> fragPile;
	/** All of the inter-fragment bonds */
	private Set<Bond> bondPile;
	/** A mapping between fragments and inter fragment bonds */
	private Map<Fragment,Set<Bond>> fragToInterFragmentBond;
	/** A builder for fragments specified as SMILES */
	private SMILESFragmentBuilder sBuilder;
	/** A builder for fragments specified as references to a CML data file */
	private CMLFragmentBuilder cmlBuilder;
	/** A source of unique integers */
	private IDManager idManager;

	/** Sets up a new Fragment mananger, containing no fragments.
	 *
	 * @param sBuilder A SMILESFragmentBuilder - dependency injection.
	 * @param cmlBuilder A CMLFragmentBuilder - dependency injection.
	 * @param idManager An IDManager.
	 *
	 * @throws Exception If the CML fragment file can't be found or otherwise used
	 */
	FragmentManager(SMILESFragmentBuilder sBuilder, CMLFragmentBuilder cmlBuilder, IDManager idManager) {
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
		Fragment newFrag = cmlBuilder.build(idStr, type, subType, idManager);
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
	 * @return The built fragment
	 * @throws StructureBuildingException
	 */
	Fragment buildSMILES(String smiles, String type, String subType, String labelMapping) throws StructureBuildingException {
		Fragment newFrag = sBuilder.build(smiles, type, subType, labelMapping, idManager);
		addFragment(newFrag);
		return newFrag;
	}

	/**Creates a new fragment, containing all of the atoms and bonds
	 * of all of the other fragments - ie the whole molecule. This is non-destructive,
	 * and does not update which fragment the Atoms think they are in. Atoms and Bonds
	 * are not copied.
	 *
	 * @return The unified fragment
	 */
	Fragment getUnifiedFragment() {
		Fragment outFrag = new Fragment();
		for(Fragment f : fragPile) {
			outFrag.importFrag(f);
		}
		for (Bond interFragmentBond : bondPile) {
			outFrag.addBond(interFragmentBond, false);//already will have been associated with atoms hence false
		}
		return outFrag;
	}

	/** Joins two fragments together, by creating a bond between them
	 *
	 * @param fromAtom The identity of an atom on one fragment
	 * @param toAtom The identity of an atom on another fragment
	 * @param bondOrder The order of the joining bond
	 */
	void attachFragments(Atom fromAtom, Atom toAtom, int bondOrder) {
		Bond b =new Bond(fromAtom, toAtom, bondOrder);
		bondPile.add(b);
		fragToInterFragmentBond.get(fromAtom.getFrag()).add(b);
		fragToInterFragmentBond.get(toAtom.getFrag()).add(b);
		fromAtom.addBond(b);
		toAtom.addBond(b);
	}

	/** Incorporates a fragment, usually a suffix, into a parent fragment, creating a bond between them.
	 *
	 * @param suffixFrag The fragment to be incorporated
	 * @param fromID An id on that fragment
	 * @param toFrag The parent fragment
	 * @param toID An id on that fragment
	 * @param bondOrder The order of the joining bond
	 */
	void incorporateFragment(Fragment suffixFrag, int fromID, Fragment toFrag, int toID, int bondOrder) throws StructureBuildingException {
		toFrag.importFrag(suffixFrag);
		toFrag.addBond(new Bond(toFrag.getAtomByIDOrThrow(fromID), toFrag.getAtomByIDOrThrow(toID), bondOrder));
		removeFragment(suffixFrag);
	}

	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtomID The id of the lower-numbered atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
	 */
	void unsaturate(int fromAtomID, int bondOrder, Fragment fragment) throws StructureBuildingException {
		int toAtomID = fromAtomID + 1;
		if (fragment.getAtomByID(toAtomID)==null || fragment.getAtomByID(toAtomID).getType().equals("suffix")){//allows something like cyclohexan-6-ene, something like butan-4-ene will still fail
			List<Atom> neighbours =fragment.getAtomByIDOrThrow(fromAtomID).getAtomNeighbours();
			if (neighbours.size() >=2){
				int firstID =fragment.getIdOfFirstAtom();
				for (Atom a : neighbours) {
					if (a.getID() ==firstID){
						toAtomID=firstID;
						break;
					}
				}
			}
		}
		Bond b = fragment.findBondOrThrow(fromAtomID, toAtomID);
		b.setOrder(bondOrder);
	}

	/** Adjusts the order of a bond in a fragment.
	 *
	 * @param fromAtomID The id of the first atom in the bond
	 * @param locantTo The locant of the other atom in the bond
	 * @param bondOrder The new bond order
	 * @param fragment The fragment
	 */
	void unsaturate(int fromAtomID, String locantTo, int bondOrder, Fragment fragment) throws StructureBuildingException {
		int toAtomID = fragment.getIDFromLocantOrThrow(locantTo);
		Bond b = fragment.findBondOrThrow(fromAtomID, toAtomID);
		b.setOrder(bondOrder);
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
			Fragment f = sBuilder.build(atomSymbol, idManager);
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

	/** Works out where to put an "one", if this is unspecified. position 2 for propanone
	 * and higher, else 1. Position 2 is assumed to be 1 higher than the atomIndice given.
	 *
	 * @param fragment The fragment
	 * @return the appropriate atom indice
	 */
	int findKetoneAtomIndice(Fragment fragment, int atomIndice) {
		if(fragment.getChainLength() < 3){
			return atomIndice;
		}
		else
			if (atomIndice +1>=fragment.getAtomList().size()){
				return 1;//this probably indicates a problem with the input name but nonetheless 1 is a better answer than an indice which isn't even in the range of the fragment
			}
			else{
				return atomIndice +1;
			}
	}

	/** Gets an atom, given an id number
	 *
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
	 *
	 * @param id The id of the atom
	 * @return The atom
	 */
	Atom getAtomByIDOrThrow(int id) throws StructureBuildingException {
		Atom a = getAtomByID(id);
		if(a == null) throw new StructureBuildingException("Couldn't get atom by id");
		return a;
	}

	/**Turns all of the spare valencies in the framents into double bonds.
	 *
	 * @throws Exception
	 */
	void convertSpareValenciesToDoubleBonds() throws StructureBuildingException {
		for(Fragment f : fragPile) {
			f.convertSpareValenciesToDoubleBonds();
		}
	}

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
	void addFragment(Fragment frag)  {
		fragPile.add(frag);
		fragToInterFragmentBond.put(frag, new HashSet<Bond>());
	}

	/**
	 * Removes a fragment from the fragPile and inter fragment bonds from the bondpile.
	 * Throws an exception if fragment wasn't present
	 * @param frag
	 * @throws StructureBuildingException
	 */
	void removeFragment(Fragment frag) throws StructureBuildingException {
		if (!fragPile.remove(frag)){
			throw new StructureBuildingException("Fragment not found in fragPile");
		}
		if (fragToInterFragmentBond.get(frag) !=null){
			Set<Bond> interFragmentBondsInvolvingFragment = fragToInterFragmentBond.get(frag);
			for (Bond bond : interFragmentBondsInvolvingFragment) {
				bondPile.remove(bond);
			}
			fragToInterFragmentBond.remove(frag);
		}
	}

	/**
	 * Changes the charge on the atom described by the given id by the given amount
	 * @param Id
	 * @param charge
	 * @param fragment
	 * @throws StructureBuildingException
	 */
	void changeCharge(int Id, int charge, Fragment fragment) throws StructureBuildingException {
		Atom atom= fragment.getAtomByIDOrThrow(Id);
		int currentCharge =atom.getCharge();
		atom.setCharge(currentCharge+=charge);
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
	 * labels the atoms using new ids from the idManager and adds to the fragManager in state
	 * @param originalFragment
	 * @return
	 * @throws StructureBuildingException
	 */
	Fragment copyAndRelabel(Fragment originalFragment) throws StructureBuildingException {
		return copyAndRelabel(originalFragment, null);
	}


	/**
	 * Creates a copy of a fragment by copying data
	 * labels the atoms using new ids from the idManager and adds to the fragManager in state
	 * @param originalFragment
	 * @param stringToAddToAllLocants: typically used to append primes to all locants, can be null
	 * @return
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
			int ID = idManager.getNextID();
			ArrayList<String> newLocants = new ArrayList<String>(atom.getLocants());
			if (stringToAddToAllLocants !=null){
				for (int i = 0; i < newLocants.size(); i++) {
					newLocants.set(i, newLocants.get(i) + stringToAddToAllLocants);
				}
			}
			Atom newAtom =new Atom(ID, newLocants, atom.getElement(), newFragment);
			newAtom.setCharge(atom.getCharge());
			newAtom.setSpareValency(atom.getSpareValency());
			if (atom.getAtomParityElement() != null){
				newAtom.setAtomParityElement(new Element(atom.getAtomParityElement()));
			}
			newAtom.setExplicitHydrogens(atom.getExplicitHydrogens());
			newAtom.setValency(atom.getValency());
			//outValency is derived from the outIDs so is automatically cloned
			newAtom.setAtomIsInACycle(atom.getAtomIsInACycle());
			newAtom.setType(atom.getType());//may be different from fragment type if the original atom was formerly in a suffix
			newAtom.setNotes(new HashMap<String, String>(atom.getNotes()));
			newFragment.addAtom(newAtom);
			idMap.put(atom.getID(),ID);
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
		List<Bond> bondList =originalFragment.getBondList();
		for (Bond bond : bondList) {
			Bond newBond=new Bond(newFragment.getAtomByIDOrThrow(idMap.get(bond.getFrom())),newFragment.getAtomByIDOrThrow(idMap.get(bond.getTo())), bond.getOrder());
			newBond.setSmilesStereochemistry(bond.getSmilesStereochemistry());
			if (bond.getBondStereoElement() != null){
				newBond.setBondStereoElement(new Element(bond.getBondStereoElement()));
			}
			newFragment.addBond(newBond);
		}
		addFragment(newFragment);
		return newFragment;
	}

	static void relabelFusedRingSystem(Fragment fusedring){
		relabelFusedRingSystem(fusedring.getAtomList());
	}

	/**Adjusts the labeling on a fused ring system, such that bridgehead atoms
	 * have locants endings in 'a' or 'b' etc. Example: naphthalene
	 * 1,2,3,4,5,6,7,8,9,10->1,2,3,4,4a,5,6,7,8,8a
	 */
	static void relabelFusedRingSystem(List<Atom> atomList) {
		int locantVal = 0;
		char locantLetter = 'a';
		for (Atom atom : atomList) {
			atom.clearLocants();
		}
		for (Atom atom : atomList) {
			if(!atom.getElement().equals("C") || atom.getBonds().size() < 3) {
				locantVal++;
				locantLetter = 'a';
				atom.addLocant(Integer.toString(locantVal));
			} else {
				atom.addLocant(Integer.toString(locantVal) + locantLetter);
				locantLetter++;
			}
		}
	}

	/**
	 * Assign element locants to groups/suffixes. These are in addition to any numerical locants that are present.
	 * Adds primes to make each locant unique.
	 * For groups a locant is not given to carbon atoms
	 * If an element appears in a suffix then element locants are not assigned to occurrences of that element in the parent group
	 * @param suffixableFragment
	 * @param suffixFragments
	 */
	static void assignElementLocants(Fragment suffixableFragment, ArrayList<Fragment> suffixFragments) {
		HashMap<String,Integer> elementCount =new HashMap<String,Integer>();//keeps track of how many times each element has been seen

		HashSet<Atom> atomsToIgnore = new HashSet<Atom>();//atoms which already have a symbolic locant
		ArrayList<Fragment> allFragments =new ArrayList<Fragment>(suffixFragments);
		allFragments.add(suffixableFragment);
		/*
		 * First check whether any element locants have already been assigned, these will take precedence
		 */
		for (Fragment fragment : allFragments) {
			List<Atom> atomList =fragment.getAtomList();
			for (Atom atom : atomList) {
				List<String> elementSymbolLocants =atom.getElementSymbolLocants();
				for (String locant : elementSymbolLocants) {
					int primeCount =0;
					for(int i=0;i<locant.length();i++) {
						if(locant.charAt(i) == '\'') primeCount++;
					}
					String element =locant.substring(0, locant.length()-primeCount);
					if (elementCount.get(element)==null || (elementCount.get(element) < primeCount +1)){
						elementCount.put(element, primeCount +1);
					}
					atomsToIgnore.add(atom);
				}
			}
		}

		for (Fragment fragment : suffixFragments) {
			List<Atom> atomList =fragment.getAtomList();

			/*
			 * Sort them such that single bonded atoms are higher priority than doubled bonded atoms
			 */
			Collections.sort(atomList, new SortAtomsForElementSymbols());
			for (Atom atom : atomList) {//add the locants
				if (atomsToIgnore.contains(atom)){continue;}
				String element =atom.getElement();
				if (elementCount.get(element)==null){
					atom.addLocant(element);
					elementCount.put(element,1);
				}
				else{
					int count =elementCount.get(element);
					atom.addLocant(element + StringTools.multiplyString("'", count));
					elementCount.put(element, count +1);
				}
			}
		}
		HashSet<String> elementToIgnore = new HashSet<String>(elementCount.keySet());
		elementCount =new HashMap<String,Integer>();
		List<Atom> atomList =suffixableFragment.getAtomList();
		Atom atomToAddCLabelTo=null;//only add a C label if there is only one C in the main group
		for (Atom atom : atomList) {
			if (atomsToIgnore.contains(atom)){continue;}
			String element =atom.getElement();
			if (elementToIgnore.contains(element)){
				continue;
			}
			if (element.equals("C")){
				if (atomToAddCLabelTo !=null){
					elementToIgnore.add("C");
					atomToAddCLabelTo=null;
				}
				else{
					atomToAddCLabelTo =atom;
				}
				continue;
			}
			if (elementCount.get(element)==null){
				atom.addLocant(element);
				elementCount.put(element,1);
			}
			else{
				int count =elementCount.get(element);
				atom.addLocant(element + StringTools.multiplyString("'", count));
				elementCount.put(element, count +1);
			}
		}
		if (atomToAddCLabelTo !=null){
			atomToAddCLabelTo.addLocant("C");
		}
	}

	/**
	 * Takes an element and produces a copy of it. Groups and suffixes are copied so that the new element
	 * has it's own group and suffix fragments
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
	 * has it's own group and suffix fragments
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
			attachFragments(fromAtom, toAtom, bond.getOrder());
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
		if (atomNeighbours.size() > 1 || (atomNeighbours.size() == 1 &&  terminalAtom.getOutValency() >0)){
			throw new StructureBuildingException("Atom to be attached to fragment should only have one bond");
		}
		terminalAtom.setElement(atomThatReplacesTerminal.getElement());
		terminalAtom.setCharge(atomThatReplacesTerminal.getCharge());
		terminalAtom.clearLocants();
		for (String locant : atomThatReplacesTerminal.getLocants()) {
			terminalAtom.addLocant(locant);
		}

		parentFrag.importFrag(childFrag);
		List<Atom> neighbours = atomThatReplacesTerminal.getAtomNeighbours();
		for (Atom neighbour : neighbours) {
			parentFrag.addBond(new Bond(terminalAtom, neighbour, childFrag.findBond(atomThatReplacesTerminal, neighbour).getOrder()));
		}
		parentFrag.removeAtom(atomThatReplacesTerminal, this);
		removeFragment(childFrag);
	}

	/**
	 * Checks if this bond is an inter fragment bond and if it is removes it
	 * @param bond
	 */
	void removeInterFragmentBondIfPresent(Bond bond) {
		if (bondPile.remove(bond)){
			fragToInterFragmentBond.get(bond.getFromAtom().getFrag()).remove(bond);
			fragToInterFragmentBond.get(bond.getToAtom().getFrag()).remove(bond);
		}
	}

	/**
	 * Gets a set of the inter fragment bonds a fragment is involved in
	 * @param frag
	 * @return set of inter fragment bonds
	 */
	Set<Bond> getInterFragmentBonds(Fragment frag) {
		return fragToInterFragmentBond.get(frag);
	}
}