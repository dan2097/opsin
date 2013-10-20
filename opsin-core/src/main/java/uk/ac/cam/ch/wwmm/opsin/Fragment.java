package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;

import nu.xom.Attribute;
import nu.xom.Element;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**A fragment of a molecule, holds bonds and atoms.
 *
 * @author ptc24
 * @author dl387
 *
 */
class Fragment {

	/**A mapping between IDs and the atoms in this fragment, by default is ordered by the order atoms are added to the fragment*/
	private final LinkedHashMap<Integer, Atom> atomMapFromId = new LinkedHashMap<Integer, Atom>();

	/**Equivalent to and synced to atomMapFromId.values() */
	private final Collection<Atom> atomCollection = atomMapFromId.values();

	/**A mapping between locants and the atoms in this fragment*/
	private final HashMap<String, Atom> atomMapFromLocant = new HashMap<String, Atom>();

	/**The bonds in the fragment*/
	private final Set<Bond> bondSet = new LinkedHashSet<Bond>();

	/**The type of the fragment, for the purpose of resolving suffixes*/
	private String type = "";

	/**The subType of the fragment, for the purpose of resolving suffixes*/
	private String subType = "";

	/**The atoms that are used when this fragment is connected to another fragment. Unused outAtoms means that the fragment is a radical or an error has occurred
	 * Initially empty */
	private final LinkedList<OutAtom> outAtoms = new LinkedList<OutAtom>();

	/**The atoms that are used on this fragment to form things like esters
	 * Initially empty */
	private final LinkedList<FunctionalAtom> functionalAtoms = new LinkedList<FunctionalAtom>();

	/**The atom that fragments connecting to this fragment connect to if a locant has not been specified
	 * Defaults to the first atom to be added to the fragment. This is typically the one with locant 1
	 * but allows for fragments with no locants. Can be overridden*/
	private Atom defaultInAtom = null;

	/**The atoms in the fragment that have been indicated to have hydrogen at the SMILES level.*/
	private final List<Atom> indicatedHydrogen = new ArrayList<Atom>();

	/**DO NOT CALL DIRECTLY EXCEPT FOR TESTING
	 * Makes an empty Fragment with a given type and subType.
	 * @param type The type of the fragment
	 * @param subType The subtype of the fragment 
	 * @throws StructureBuildingException
	 */
	Fragment(String type, String subType) throws StructureBuildingException {
		if (type==null){
			throw new StructureBuildingException("Type specified for fragment is null");
		}
		if (subType==null){
			throw new StructureBuildingException("subType specified for fragment is null");
		}
		this.type = type;
		this.subType = subType;
	}

	/**DO NOT CALL DIRECTLY EXCEPT FOR TESTING
	 * Makes an empty fragment with no specified type.*/
	Fragment() {}

	/**Produces a CML element, corresponding to the molecule. The cml element contains
	 * a molecule, which contains an atomArray and bondArray filled with atoms and bonds.
	 * The molecule element has a dummy id of m1.
	 * @param chemicalName
	 *
	 * @return The CML element.
	 * @see Atom
	 * @see Bond
	 */
	Element toCMLMolecule(String chemicalName) {
		Element cml = new Element("cml", CML_NAMESPACE);
		cml.addAttribute(new Attribute("convention","conventions:molecular"));
		cml.addNamespaceDeclaration("conventions", "http://www.xml-cml.org/convention/");
		cml.addNamespaceDeclaration("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
		cml.addNamespaceDeclaration("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
		Element molecule = new Element("molecule", CML_NAMESPACE);
		Element name = new Element("name", CML_NAMESPACE);
		name.appendChild(chemicalName);
		name.addAttribute(new Attribute("dictRef","nameDict:unknown"));
		molecule.appendChild(name);
		molecule.addAttribute(new Attribute("id", "m1"));
		Element atomArray = new Element("atomArray", CML_NAMESPACE);
		for(Atom atom : atomCollection) {
			atomArray.appendChild(atom.toCMLAtom());
		}
		Element bondArray = new Element("bondArray", CML_NAMESPACE);
		for(Bond bond : bondSet) {
			bondArray.appendChild(bond.toCMLBond());
		}
		molecule.appendChild(atomArray);
		molecule.appendChild(bondArray);
		cml.appendChild(molecule);
		return cml;
	}

	/**Adds an atom to the fragment and associates it with this fragment*/
	void addAtom(Atom atom) {
		if (defaultInAtom == null){//the first atom added becomes the defaultInAtom
			defaultInAtom = atom;
		}
		List<String> locants =atom.getLocants();
		for (String locant: locants) {
			atomMapFromLocant.put(locant, atom);
		}
		atomMapFromId.put(atom.getID(), atom);
		atom.setFrag(this);
	}
	
	/**
	 * Return the number of atoms in the fragment
	 * @return
	 */
	int getAtomCount() {
		return atomCollection.size();
	}

	/**
	 * Returns a copy of the fragment's atoms
	 * @return
	 */
	List<Atom> getAtomList() {
		return new ArrayList<Atom>(atomCollection);
	}

	/**
	 * Adds a bond to the fragment.
	 * @param bond
	 */
	void addBond(Bond bond) {
		bondSet.add(bond);
	}
	
	/**Removes a bond to the fragment if it is present.
    * @param bond
    * @return*/
	boolean removeBond(Bond bond) {
		return bondSet.remove(bond);
	}

	/**Gets bondSet.*/
	Set<Bond> getBondSet() {
		return Collections.unmodifiableSet(bondSet);
	}

	/**Gets the id of the atom in the fragment with the specified locant.
	 *
	 * @param locant The locant to look for
	 * @return The id of the found atom, or 0 if it is not found
	 * @throws StructureBuildingException 
	 */
	int getIDFromLocant(String locant) throws StructureBuildingException {
		Atom a = getAtomByLocant(locant);
		if (a != null){
			return a.getID();
		}
		return 0;
	}

	/**Gets the id of the atom in the fragment with the specified locant, throwing if this fails.
	 *
	 * @param locant The locant to look for
	 * @return The id of the found atom
     * @throws StructureBuildingException
	 */
	int getIDFromLocantOrThrow(String locant) throws StructureBuildingException {
		int id = getIDFromLocant(locant);
		if(id == 0) {
			throw new StructureBuildingException("Couldn't find id from locant " + locant + ".");
		}
		return id;
	}

	/**Gets the atom in the fragment with the specified locant.
	 *
	 * @param locant The locant to look for
	 * @return The found atom, or null if it is not found
	 * @throws StructureBuildingException 
	 */
	Atom getAtomByLocant(String locant) throws StructureBuildingException {
		Atom a =atomMapFromLocant.get(locant);
		if (a != null){
			return a;
		}
		Matcher m =MATCH_AMINOACID_STYLE_LOCANT.matcher(locant);
		if (m.matches()){//e.g. N5
			Atom backboneAtom =atomMapFromLocant.get(m.group(3));//the atom corresponding to the numeric or greek component
			if (backboneAtom==null){
				return null;
			}
			a = FragmentTools.getAtomByAminoAcidStyleLocant(backboneAtom,  m.group(1), m.group(2));
			if (a != null){
				return a;
			}
		}
		return null;
	}

	/**Gets the atom in the fragment with the specified locant, throwing if this fails.
	 *
	 * @param locant The locant to look for
	 * @return The found atom
	 * @throws StructureBuildingException
	 */
	Atom getAtomByLocantOrThrow(String locant) throws StructureBuildingException {
		Atom a = getAtomByLocant(locant);
		if(a == null) {
			throw new StructureBuildingException("Could not find the atom with locant " + locant + ".");
		}
		return a;
	}

	/**Gets the atom in the fragment with the specified ID.
	 *
	 * @param id The id of the atom.
	 * @return The found atom, or null.
	 */
	Atom getAtomByID(int id) {
		return atomMapFromId.get(id);
	}

	/**Gets the atom in the fragment with the specified ID, throwing if this fails.
	 *
	 * @param id The id of the atom.
	 * @return The found atom
     * @throws StructureBuildingException
	 */
	Atom getAtomByIDOrThrow(int id) throws StructureBuildingException {
		Atom a = getAtomByID(id);
		if(a == null) {
			throw new StructureBuildingException("Couldn't find atom with id " + id + ".");
		}
		return a;
	}

	/**Finds a bond between two specified atoms the first of which must be within the fragment
	 *
	 * @param ID1 The id of one atom
	 * @param ID2 The id of the other atom
	 * @return The bond found, or null
	 */
	Bond findBond(int ID1, int ID2) {
		Atom a = atomMapFromId.get(ID1);
		if (a!=null){
			for (Bond b : a.getBonds()) {
				if((b.getFrom() == ID1 && b.getTo() == ID2) ||
						(b.getTo() == ID1 && b.getFrom() == ID2)) {
					return b;
				}
			}
		}
		return null;
	}

	/**Finds a bond between two specified atoms the first of which must be within the fragment, throwing if it fails.
	 *
	 * @param ID1 The id of one atom
	 * @param ID2 The id of the other atom
	 * @return The bond found
     * @throws StructureBuildingException
	 */
	Bond findBondOrThrow(int ID1, int ID2) throws StructureBuildingException {
		Bond b = findBond(ID1, ID2);
		if(b == null) {
			throw new StructureBuildingException("Couldn't find specified bond");
		}
		return b;
	}

	/**Works out how many atoms there are in the fragment there are
	 * with consecutive locants, starting from 1 that are in a chain
	 *
	 * @return The number of atoms in the locant chain.
	 * @throws StructureBuildingException 
	 */
	int getChainLength() throws StructureBuildingException {
		int length = 0;
		Atom next = getAtomByLocant(Integer.toString(length+1));
		Atom previous = null;
		while (next !=null){
			if (previous !=null && previous.getBondToAtom(next) == null){
				break;
			}
			length++;
			previous = next;
			next = getAtomByLocant(Integer.toString(length+1));
		}
		return length;
	}

	/**Gets type.*/
	String getType() {
		return type;
	}

	/**Gets subType.
    * @return subType*/
	String getSubType() {
		return subType;
	}

	/**
	 * How many OutAtoms (i.e. radicals) are associated with this fragment
	 * @return
	 */
	int getOutAtomCount() {
		return outAtoms.size();
	}

	/**
	 * Gets the outAtom at a specific index of the outAtoms linkedList
	 * @param i
	 * @return
	 */
	OutAtom getOutAtom(int i) {
		return outAtoms.get(i);
	}

	/**
	 * Adds an outAtom
	 * @param id
	 * @param valency
	 * @param setExplicitly
	 * @throws StructureBuildingException 
	 */
	void addOutAtom(int id, int valency, Boolean setExplicitly) throws StructureBuildingException {
		addOutAtom(getAtomByIDOrThrow(id), valency, setExplicitly);
	}
	
	/**
	 * Adds an outAtom
	 * @param atom
	 * @param valency
	 * @param setExplicitly
	 */
	void addOutAtom(Atom atom, int valency, Boolean setExplicitly) {
		outAtoms.add(new OutAtom(atom, valency, setExplicitly));
	}
	
	/**
	 * Includes the OutAtoms of a given fragment into this fragment
	 * Note that no OutAtoms are created in doing this 
	 * @param outAtoms
	 */
	void incorporateOutAtoms(Fragment frag) {
		outAtoms.addAll(frag.outAtoms);
	}

	/**
	 * Removes the outAtom at a specific index of the outAtom linkedList
	 * @param i
	 */
	void removeOutAtom(int i) {
		OutAtom removedOutAtom = outAtoms.remove(i);
		if (removedOutAtom.isSetExplicitly()){
			removedOutAtom.getAtom().addOutValency(-removedOutAtom.getValency());
		}
	}

	/**
	 * Removes the specified outAtom from the outAtoms linkedList
	 * @param outAtom
	 */
	void removeOutAtom(OutAtom outAtom) {
		if (outAtoms.remove(outAtom) && outAtom.isSetExplicitly()){
			outAtom.getAtom().addOutValency(-outAtom.getValency());
		}
	}

	/**
	 * How many functionalAtoms (i.e. locations that can form esters) are associated with this fragment
	 * @return
	 */
	int getFunctionalAtomCount() {
		return functionalAtoms.size();
	}

	/**
	 * Gets the functionalAtom at a specific index of the functionalAtoms linkedList
	 * @param i
	 * @return
	 */
	FunctionalAtom getFunctionalAtom(int i) {
		return functionalAtoms.get(i);
	}
	
	/**Adds a functionalAtom
    * @param atom*/
	void addFunctionalAtom(Atom atom) {
		functionalAtoms.add(new FunctionalAtom(atom));
	}
	
	/**
	 * Includes the FunctionalAtoms of a given fragment into this fragment
	 * Note that no FunctionalAtoms are created in doing this 
	 * @param functionalAtoms
	 */
	void incorporateFunctionalAtoms(Fragment frag) {
		functionalAtoms.addAll(frag.functionalAtoms);
	}

	/**
	 * Removes the functionalAtom at a specific index of the functionalAtoms linkedList
	 * @param i
	 * @return 
	 */
	FunctionalAtom removeFunctionalAtom(int i) {
		return functionalAtoms.remove(i);
	}

	/**
	 * Removes the specified functionalAtom from the functionalAtoms linkedList
	 * @param functionalAtom
	 */
	void removeFunctionalAtom(FunctionalAtom functionalAtom) {
		functionalAtoms.remove(functionalAtom);
	}

	/**Gets a list of atoms in the fragment that connect to a specified atom
	 *
	 * @param atom The reference atom
	 * @return The list of atoms connected to the atom
	 * @throws StructureBuildingException
	 */
	List<Atom> getIntraFragmentAtomNeighbours(Atom atom) throws StructureBuildingException {
		List<Atom> results = new ArrayList<Atom>();
		for(Bond b :  atom.getBonds()) {
			//recalled atoms will be null if they are not part of this fragment
			if(b.getFromAtom() == atom) {
				Atom a =getAtomByID(b.getTo());
				if (a!=null){
					results.add(a);
				}
			} else if(b.getToAtom() == atom) {
				Atom a =getAtomByID(b.getFrom());
				if (a!=null){
					results.add(a);
				}
			}
			else{
				throw new StructureBuildingException("A bond associated with an atom does not involve it");
			}
		}
		return results;
	}

	/**Calculates the number of bonds connecting to the atom, excluding bonds to implicit
	 * hydrogens. Double bonds count as
	 * two bonds, etc. Eg ethene - both C's have an incoming valency of 2.
	 *
	 * Only bonds to atoms within the fragment are counted. Suffix atoms are excluded
	 *
	 * @param atom
     * @return Incoming Valency
	 * @throws StructureBuildingException
	 */
	int getIntraFragmentIncomingValency(Atom atom) throws StructureBuildingException {
		int v = 0;
		for(Bond b :  atom.getBonds()) {
			//recalled atoms will be null if they are not part of this fragment
			if(b.getFromAtom() == atom) {
				Atom a =getAtomByID(b.getTo());
				if (a!=null && !a.getType().equals(SUFFIX_TYPE_VAL)){
					v += b.getOrder();
				}
			} else if(b.getToAtom() == atom) {
				Atom a =getAtomByID(b.getFrom());
				if (a!=null && !a.getType().equals(SUFFIX_TYPE_VAL)){
					v += b.getOrder();
				}
			}
			else{
				throw new StructureBuildingException("A bond associated with an atom does not involve it");
			}
		}
		return v;
	}

	/**
	 * Checks valencies are sensible
	 * @throws StructureBuildingException
	 */
	void checkValencies() throws StructureBuildingException {
		for(Atom a : atomCollection) {
			if(!ValencyChecker.checkValency(a)) {
				throw new StructureBuildingException("Atom is in unphysical valency state! Element: " + a.getElement() + " valency: " + a.getIncomingValency());
			}
		}
	}

	/**
	 * Removes an atom from this fragment
	 * @param atom
	 */
	void removeAtom(Atom atom) {
		int atomID =atom.getID();
		atomMapFromId.remove(atomID);
		for (String l : atom.getLocants()) {
			atomMapFromLocant.remove(l);
		}
		if (defaultInAtom == atom){
			defaultInAtom = getFirstAtom();
		}
	}
	/**
	 * Retrieves the overall charge of the fragment by querying all its atoms
	 * @return
	 */
	int getCharge() {
		int charge=0;
		for (Atom a : atomCollection) {
			charge+=a.getCharge();
		}
		return charge;
	}

	/**
	 * Sets the type of the fragment e.g. aromaticStem
	 * @param type
	 */
	void setType(String type) {
		this.type = type;
	}
	
	/**
	 * Sets the subType of the fragment
	 * @param subType
	 */
	void setSubType(String subType) {
		this.subType = subType;
	}

	Atom getDefaultInAtom() {
		return defaultInAtom;
	}

	void setDefaultInAtom(Atom inAtom) {
		this.defaultInAtom=inAtom;
	}

	/**
	 * Adds a mapping between the locant and atom object
	 * @param locant A locant as a string
	 * @param a An atom
	 */
	void addMappingToAtomLocantMap(String locant, Atom a){
		atomMapFromLocant.put(locant, a);
	}

	/**
	 * Removes a mapping between a locant
	 * @param locant A locant as a string
	 */
	void removeMappingFromAtomLocantMap(String locant){
		atomMapFromLocant.remove(locant);
	}

	/**
	 * Checks to see whether a locant is present on this fragment
	 * @param locant
	 * @return
	 * @throws StructureBuildingException
	 */
	boolean hasLocant(String locant) throws StructureBuildingException {
		return getAtomByLocant(locant)!=null;
	}
	

	/**
	 * Returns an unmodifiable list of the locants associated with this fragment
	 * @return
	 */
	Set<String> getLocants() {
		return Collections.unmodifiableSet(atomMapFromLocant.keySet());
	}

	List<Atom> getIndicatedHydrogen() {
		return indicatedHydrogen;
	}
	
	void addIndicatedHydrogen(Atom atom) {
		indicatedHydrogen.add(atom);
	}


	Atom getAtomOrNextSuitableAtomOrThrow(Atom startingAtom, int additionalValencyRequired, boolean takeIntoAccountOutValency) throws StructureBuildingException {
		Atom a =getAtomOrNextSuitableAtom(startingAtom, additionalValencyRequired, takeIntoAccountOutValency);
		if (a==null){
			throw new StructureBuildingException("No suitable atom found");
		}
		return a;
	}

	/**
	 * Takes an id and additional valency required. Returns the atom associated with that id if adding the specified valency will not violate
	 * that atom type's maximum valency.
	 * If this is not possible it iterates sequentially through all atoms in the fragment till one is found
	 * Spare valency is initially taken into account so that the atom is not dearomatised
	 * If this is impossible to accomplish dearomatisation is done
	 * If an atom is still not found an exception is thrown
	 * atoms belonging to suffixes are never selected unless the original id specified was a suffix atom
	 * @param startingAtom
	 * @param additionalValencyRequired The increase in valency that will be required on the desired atom
	 * @param takeIntoAccountOutValency
     * @return Atom
	 */
	Atom getAtomOrNextSuitableAtom(Atom startingAtom, int additionalValencyRequired, boolean takeIntoAccountOutValency) {
		List<Atom> atomList =getAtomList();
		Atom currentAtom = startingAtom;
		int atomCounter=0;
		int atomListPosition=atomList.indexOf(currentAtom);
		int startingIndex =atomListPosition;

		do {//aromaticity preserved and standard valency assumed
			atomCounter++;
			if (atomListPosition >= atomList.size()){
				atomListPosition -=(atomList.size());
			}
			currentAtom=atomList.get(atomListPosition);
			if (FragmentTools.isCharacteristicAtom(currentAtom)){
				atomListPosition++;
				continue;
			}
			int currentExpectedValency = currentAtom.determineValency(takeIntoAccountOutValency);
			if(currentExpectedValency >= (currentAtom.getIncomingValency() + additionalValencyRequired + (currentAtom.hasSpareValency() ? 1 : 0) + (takeIntoAccountOutValency ? currentAtom.getOutValency() : 0))){
				return currentAtom;
			}
			atomListPosition++;
		}
		while(atomCounter < atomList.size());

		atomListPosition =startingIndex;
		atomCounter=0;

		do {//aromaticity preserved, standard valency assumed, non functional suffixes substitutable
			atomCounter++;
			if (atomListPosition >= atomList.size()){
				atomListPosition -=(atomList.size());
			}
			currentAtom=atomList.get(atomListPosition);
			if (FragmentTools.isFunctionalAtomOrAldehyde(currentAtom)){
				atomListPosition++;
				continue;
			}
			int currentExpectedValency = currentAtom.determineValency(takeIntoAccountOutValency);
			if(currentExpectedValency >= (currentAtom.getIncomingValency() + additionalValencyRequired + (currentAtom.hasSpareValency() ? 1 : 0) + (takeIntoAccountOutValency ? currentAtom.getOutValency() : 0))){
				return currentAtom;
			}
			atomListPosition++;
		}
		while(atomCounter < atomList.size());
		
		atomListPosition =startingIndex;
		atomCounter=0;

		do {//aromaticity preserved any suffix substitutable
			atomCounter++;
			if (atomListPosition >= atomList.size()){
				atomListPosition -=(atomList.size());
			}
			currentAtom=atomList.get(atomListPosition);

			if(ValencyChecker.checkValencyAvailableForBond(currentAtom, additionalValencyRequired + (currentAtom.hasSpareValency() ? 1 : 0) + (takeIntoAccountOutValency ? currentAtom.getOutValency() : 0))){
				return currentAtom;
			}
			atomListPosition++;
		}
		while(atomCounter < atomList.size());

		atomListPosition =startingIndex;
		atomCounter=0;
		do {//aromaticity dropped, anything substitutable
			atomCounter++;
			if (atomListPosition >= atomList.size()){
				atomListPosition -=(atomList.size());
			}
			currentAtom=atomList.get(atomListPosition);
			if(ValencyChecker.checkValencyAvailableForBond(currentAtom, additionalValencyRequired + (takeIntoAccountOutValency ? currentAtom.getOutValency() : 0))){
				return currentAtom;
			}
			atomListPosition++;
		}
		while(atomCounter < atomList.size());
		return null;
	}

	/**
	 * Returns the id of the first atom in the fragment
	 * @return
	 * @throws StructureBuildingException
	 */
	int getIdOfFirstAtom() {
		return getFirstAtom().getID();
	}

	/**
	 * Returns the the first atom in the fragment or null if it has no atoms
	 * Typically the first atom will be the first atom that was added to the fragment
	 * @return firstAtom
	 */
	Atom getFirstAtom(){
		Iterator<Atom> atomIterator =atomCollection.iterator();
		if (atomIterator.hasNext()){
			return atomIterator.next();
		}
		return null;
	}

	/**
	 * Clears and recreates atomMapFromId (and hence AtomCollection) using the order of the atoms in atomList
	 * @param atomList
	 * @throws StructureBuildingException
	 */
	void reorderAtomCollection(List<Atom> atomList) throws StructureBuildingException {
		if (atomMapFromId.size()!=atomList.size()){
			throw new StructureBuildingException("atom list is not the same size as the number of atoms in the fragment");
		}
		atomMapFromId.clear();
		for (Atom atom : atomList) {
			atomMapFromId.put(atom.getID(), atom);
		}
	}

	/**
	 * Reorders the fragment's internal atomList by the value of the first locant of the atoms
	 * e.g. 1,2,3,3a,3b,4
	 * Used for assuring the correct order of atom iteration when performing ring fusion
	 * @throws StructureBuildingException 
	 */
	void sortAtomListByLocant() throws StructureBuildingException {
		List<Atom> atomList =getAtomList();
		Collections.sort(atomList, new FragmentTools.SortByLocants());
		reorderAtomCollection(atomList);
	}
}



