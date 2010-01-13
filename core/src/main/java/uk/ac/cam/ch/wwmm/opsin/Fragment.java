package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Attribute;
import nu.xom.Element;

/**A fragment of a molecule, holds bonds and atoms.
 *
 * @author ptc24/dl387
 *
 */
class Fragment {

	/**A mapping between IDs and the atoms in this fragment, by default is ordered by the order atoms are added to the fragment*/
	private LinkedHashMap<Integer, Atom> atomMapFromId = new LinkedHashMap<Integer, Atom>();

	/**Equivalent to and synced to atomMapFromId.values() */
	private Collection<Atom> atomCollection = atomMapFromId.values();

	/**A mapping between locants and the atoms in this fragment*/
	private HashMap<String, Atom> atomMapFromLocant = new HashMap<String, Atom>();

	/**The bonds in the fragment*/
	private Set<Bond> bondSet = new LinkedHashSet<Bond>();

	/**The type of the fragment, for the purpose of resolving suffixes*/
	private String type = "";

	/**The subType of the fragment, for the purpose of resolving suffixes*/
	private String subType = "";

	/**The IDs of atoms that are used when this fragment is connected to another fragment. Unused outIDs means that the name is a radical or an error has occured
	 * Initially empty */
	private LinkedList<OutID> outIDs = new LinkedList<OutID>();

	/**The IDs of atoms that are used on this fragment to form things like esters
	 * Initially empty */
	private LinkedList<FunctionalID> functionalIDs = new LinkedList<FunctionalID>();

	/**The IDs of atoms that must be bonded to and the order of the bonds. Currently used in suffixes and in multiplicative nomenclature
	 * Initially empty */
	private LinkedList<InID> inIDs = new LinkedList<InID>();

	/**The ID of the atom that fragments connecting to the fragment connect to if a locant has not been specified
	 * Set by the first atom to be added to the fragment. This is typically the one with locant 1
	 * but allows for fragments with no locants. Can be overridden*/
	private int defaultInID = 0;

	/**The id of the indicated Hydrogen on the fragment.
	 * Defaults to null.*/
	private Integer indicatedHydrogen;

	private static Pattern matchAminoAcidStyleLocant =Pattern.compile("([A-Z][a-z]?)('*)(\\d+[a-z]?'*)");
	private static Pattern matchNumericLocant =Pattern.compile("\\d+[a-z]?'*");

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

	/**Produces a CML cml element, corresponding to the molecule. The cml element contains
	 * a molecule, which contains an atomArray and bondArray filled with atoms and bonds.
	 * The molecule element has a dummy id of m1.
	 * @param chemicalName
	 *
	 * @return The CML element.
	 * @see Atom
	 * @see Bond
	 */
	Element toCMLMolecule(String chemicalName) {
		Element cml = new Element("cml");
		cml.addAttribute(new Attribute("convention","cmlDict:cmllite"));
		cml.addNamespaceDeclaration("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
		cml.addNamespaceDeclaration("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
		Element molecule = new Element("molecule");
		Element name = new Element("name");
		name.appendChild(chemicalName);
		name.addAttribute(new Attribute("dictRef","nameDict:unknown"));
		molecule.appendChild(name);
		molecule.addAttribute(new Attribute("id", "m1"));
		Element atomArray = new Element("atomArray");
		for(Atom atom : atomCollection) {
			atomArray.appendChild(atom.toCMLAtom());
		}
		Element bondArray = new Element("bondArray");
		for(Bond bond : bondSet) {
			bondArray.appendChild(bond.toCMLBond());
		}
		molecule.appendChild(atomArray);
		molecule.appendChild(bondArray);
		cml.appendChild(molecule);
		XOMTools.setNamespaceURIRecursively(cml, "http://www.xml-cml.org/schema");
		return cml;
	}

	/**Adds an atom to the fragment and associates it with this fragment*/
	void addAtom(Atom atom) {
		if (defaultInID==0){//the first atom added becomes the defaultInID
			defaultInID =atom.getID();
		}
		List<String> locants =atom.getLocants();
		for (String locant: locants) {
			atomMapFromLocant.put(locant, atom);
		}
		atomMapFromId.put(atom.getID(), atom);
		atom.setFrag(this);
	}

	/**Gets atomList.*/
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
	
	/**Removes a bond to the fragment if it is present.*/
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
	 */
	int getIDFromLocant(String locant) {
		Atom a =atomMapFromLocant.get(locant);
		if (a != null){
			return a.getID();
		}
		return 0;
	}

	/**Gets the id of the atom in the fragment with the specified locant, throwing if this fails.
	 *
	 * @param locant The locant to look for
	 * @return The id of the found atom
	 */
	int getIDFromLocantOrThrow(String locant) throws StructureBuildingException {
		int id = getIDFromLocant(locant);
		if(id == 0) throw new StructureBuildingException("Couldn't find id from locant " + locant + ".");
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
		Matcher m =matchAminoAcidStyleLocant.matcher(locant);
		if (m.matches()){//e.g. N5
			a =atomMapFromLocant.get(m.group(3));//the atom corresponding to the numeric component
			if (a==null){
				return null;
			}
			String elementSymbol = m.group(1);
			locant = elementSymbol +m.group(2);//the element symbol +primes
			//depth first search for the atom with this locant terminating whenever a numeric locant is encountered
			LinkedList<Atom> stack = new LinkedList<Atom>();
			stack.add(a);
			Set<Atom> atomsVisited =new HashSet<Atom>();
			Atom potentialMatch = null;//if you are looking for N' and you find N if no N can be found the N' will be returned
			while (stack.size() > 0) {
				Atom currentAtom =stack.removeLast();
				atomsVisited.add(currentAtom);
				List<Atom> neighbours = currentAtom.getAtomNeighbours();
				mainLoop: for (Atom neighbour : neighbours) {
					if (atomsVisited.contains(neighbour)){//already visited
						continue;
					}
					List<String> locants = neighbour.getLocants();
					if (!neighbour.getType().equals("suffix")){
						for (String neighbourLocant : locants) {
							if (matchNumericLocant.matcher(neighbourLocant).matches()){//gone to an inappropriate atom
								continue mainLoop;
							}
						}
					}
					for (String neighbourLocant : locants) {
						if (neighbourLocant.startsWith(elementSymbol)){
							if (neighbourLocant.equals(locant)){
								return neighbour;
							}
							else{
								potentialMatch =neighbour;
							}
						}
					}
					stack.add(neighbour);
				}
			}
			if (potentialMatch != null){
				return potentialMatch;
			}
		}
		return null;
	}

	/**Gets the atom in the fragment with the specified locant, throwing if this fails.
	 *
	 * @param locant The locant to look for
	 * @return The found atom
	 */
	Atom getAtomByLocantOrThrow(String locant) throws StructureBuildingException {
		Atom a = getAtomByLocant(locant);
		if(a == null) throw new StructureBuildingException("Could not find the atom with locant " + locant + ".");
		return a;
	}

	/**Gets the atom in the fragment with the specified ID.
	 *
	 * @param id The id of the atom.
	 * @return The found atom, or null.
	 */
	Atom getAtomByID(int id) {
		Atom a = atomMapFromId.get(id);
		if (a !=null){
			return a;
		}
		return null;
	}

	/**Gets the atom in the fragment with the specified ID, throwing if this fails.
	 *
	 * @param id The id of the atom.
	 * @return The found atom
	 */
	Atom getAtomByIDOrThrow(int id) throws StructureBuildingException {
		Atom a = getAtomByID(id);
		if(a == null) throw new StructureBuildingException("Couldn't find atom with id " + id + ".");
		return a;
	}

	/**Finds a bond between two specified atoms within the fragment.
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

	/**Finds a bond between two specified atoms within the fragment, throwing if fails.
	 *
	 * @param ID1 The id of one atom
	 * @param ID2 The id of the other atom
	 * @return The bond found
	 */
	Bond findBondOrThrow(int ID1, int ID2) throws StructureBuildingException {
		Bond b = findBond(ID1, ID2);
		if(b == null) throw new StructureBuildingException("Couldn't find specified bond");
		return b;
	}

	/**Finds a bond between two specified atoms within the fragment.
	 *
	 * @param a1 An atom
	 * @param a2 Another atom
	 * @return The bond found, or null
	 */
	Bond findBond(Atom a1, Atom a2) {
		for (Bond b : a1.getBonds()) {
			if((b.getFromAtom() == a1 && b.getToAtom() == a2) ||
					(b.getToAtom() == a1 && b.getFromAtom() == a2)) {
				return b;
			}
		}
		return null;
	}

	/**Finds a bond between two specified atoms within the fragment, throwing if fails.
	 *
	 * @param a1 An atom
	 * @param a2 Another atom
	 * @return The bond found
	 * @throws StructureBuildingException
	 */
	Bond findBondOrThrow(Atom a1, Atom a2) throws StructureBuildingException {
		Bond b = findBond(a1, a2);
		if(b == null) throw new StructureBuildingException("Couldn't find specified bond");
		return b;
	}

	/**Works out how many atoms there are in the fragment there are
	 * with consecutive locants, starting from 1;
	 *
	 * @return The number of atoms in the locant chain.
	 */
	int getChainLength() {
		int length=0;
		while(getIDFromLocant(Integer.toString(length+1)) > 0) {
			length++;
		}
		return length;
	}

	/**Gets type.*/
	String getType() {
		return type;
	}

	/**Gets subType.*/
	String getSubType() {
		return subType;
	}

	/**Gets the linkedList of outIDs. This is not modifiable, use the relevant methods in this class to modify it*/
	List<OutID> getOutIDs() {
		return Collections.unmodifiableList(outIDs);
	}

	/**
	 * Gets the outID at a specific index of the outIDs linkedList
	 * @param i
	 * @return
	 */
	OutID getOutID(int i) {
		return outIDs.get(i);
	}

	/**Adds an outID
	 * @throws StructureBuildingException */
	void addOutID(int id, int valency, Boolean setExplicitly) throws StructureBuildingException {
		outIDs.add(new OutID(id, valency, setExplicitly, this));
		if (setExplicitly){
			getAtomByIDOrThrow(id).addOutValency(valency);
		}
	}

	/**Adds a list of outIDs, copies of the given outIDs are not made*/
	void addOutIDs(List<OutID> outIDs) {
		this.outIDs.addAll(outIDs);
	}

	/**
	 * Removes the outID at a specific index of the outIDs linkedList
	 * @param i
	 * @throws StructureBuildingException
	 */
	void removeOutID(int i) throws StructureBuildingException {
		OutID removedOutID = outIDs.remove(i);
		if (removedOutID.setExplicitly){
			getAtomByIDOrThrow(removedOutID.id).addOutValency(-removedOutID.valency);
		}
	}

	/**
	 * Removes the specified outID from the outIDs linkedList
	 * @param outID
	 * @throws StructureBuildingException
	 */
	void removeOutID(OutID outID) throws StructureBuildingException {
		outIDs.remove(outID);
		if (outID.setExplicitly){
			getAtomByIDOrThrow(outID.id).addOutValency(-outID.valency);
		}
	}

	/**Gets the linkedList of inIDs. This is not modifiable, use the relevant methods in this class to modify it*/
	List<InID> getInIDs() {
		return Collections.unmodifiableList(inIDs);
	}

	/**
	 * Gets the inID at a specific index of the inIDs linkedList
	 * @param i
	 * @return
	 */
	InID getInID(int i) {
		return inIDs.get(i);
	}

	/**Adds an inID
	 * @throws StructureBuildingException */
	void addInID(int id, int valency) throws StructureBuildingException {
		inIDs.add(new InID(id, valency, this));
		getAtomByIDOrThrow(id).addOutValency(valency);
	}
	
	/**Adds a list of inIDs, copies of the given inIDs are not made*/
	void addInIDs(List<InID> inIDs) {
		this.inIDs.addAll(inIDs);
	}

	/**
	 * Removes the inID at a specific index of the inIDs linkedList
	 * @param i
	 * @throws StructureBuildingException
	 */
	void removeInID(int i) throws StructureBuildingException {
		InID removedinID = inIDs.remove(i);
		getAtomByIDOrThrow(removedinID.id).addOutValency(-removedinID.valency);
	}

	/**
	 * Removes the specified inID from the inIDs linkedList
	 * @param inID
	 * @throws StructureBuildingException
	 */
	void removeInID(InID inID) throws StructureBuildingException {
		inIDs.remove(inID);
		getAtomByIDOrThrow(inID.id).addOutValency(-inID.valency);
	}

	/**Gets the linkedList of functionalIDs*/
	List<FunctionalID> getFunctionalIDs() {
		return Collections.unmodifiableList(functionalIDs);
	}

	/**
	 * Gets the functionalID at a specific index of the functionalIDs linkedList
	 * @param i
	 * @return
	 */
	FunctionalID getFunctionalID(int i) {
		return functionalIDs.get(i);
	}

	/**Adds a functionalID*/
	void addFunctionalID(int id) {
		functionalIDs.add(new FunctionalID(id, this));
	}

	/**Adds a list of functionalIDs, copies of the given functionalIDs are not made*/
	void addFunctionalIDs(List<FunctionalID> functionalIDs) {
		this.functionalIDs.addAll(functionalIDs);
	}

	/**
	 * Removes the functionalID at a specific index of the functionalIDs linkedList
	 * @param i
	 */
	void removeFunctionalID(int i) {
		functionalIDs.remove(i);
	}

	/**
	 * Removes the specified functionalID from the functionalIDs linkedList
	 * @param functionalID
	 */
	void removeFunctionalID(FunctionalID functionalID) {
		functionalIDs.remove(functionalID);
	}

	/** Looks for the atom that would have had a hydrogen indicated,
	 * adds a spareValency to that atom, and sets indicatedHydrogen.
	 */
	void pickUpIndicatedHydrogen() {
		int svCount = 0;
		int dvCount = 0; /* Divalent, like O */
		for(Atom a : atomCollection) {
			svCount += a.hasSpareValency() ? 1 : 0;
			if(a.getElement().equals("O")) dvCount++;
			if(a.getElement().equals("S")) dvCount++;
			if(a.getElement().equals("Se")) dvCount++;
			if(a.getElement().equals("Te")) dvCount++;
		}
		if(svCount == atomCollection.size() - (1 + dvCount) && svCount > 0) {
			// Now we're looking for a trivalent C, or failing that a divalent N/P/As/Sb
			Atom nCandidate = null;
			for(Atom a : atomCollection) {
				if(a.getAtomIsInACycle() && !a.hasSpareValency()) {
					String element =a.getElement();
					if (element.equals("C")){
						indicatedHydrogen = a.getID();
						a.setSpareValency(true);
						return;
					} else if(element.equals("N") || element.equals("P") || element.equals("As") || element.equals("Sb")) {
						nCandidate = a;
					}
				}
			}
			if(nCandidate != null) {
				indicatedHydrogen = nCandidate.getID();
				nCandidate.setSpareValency(true);
			}
		}
	}

	/**Looks for double and higher bonds, converts them to single bonds
	 * and adds corresponding spareValencies to the atoms they join.
	 */
	void convertHighOrderBondsToSpareValencies()  {
		for(Bond b : bondSet) {
			if(b.getOrder() == 2) {
				Atom firstAtom =b.getFromAtom();
				Atom secondAtom =b.getToAtom();
				if (firstAtom.getAtomIsInACycle() && secondAtom.getAtomIsInACycle()){
					b.setOrder(1);
					firstAtom.setSpareValency(true);
					secondAtom.setSpareValency(true);
				}
			}
		}
		pickUpIndicatedHydrogen();
	}

	/**Increases the order of bonds joining atoms with spareValencies,
	 * and uses up said spareValencies.
	 * @throws StructureBuildingException If the algorithm can't work out where to put the bonds
	 */
	void convertSpareValenciesToDoubleBonds() throws StructureBuildingException {
		/* pick atom, getAtomNeighbours, decideIfTerminal, resolve */

		/*
		 * Correct spare valency by looking at valencyState of atom
		 *
		 */
		for(Atom a : atomCollection) {
			a.ensureSVIsConsistantWithValency(true);
		}

		/*
		 * Remove spare valency on atoms which may not form higher order bonds
		 */
		atomLoop: for(Atom a : atomCollection) {
			if(a.hasSpareValency()) {
				for(Atom aa : getAtomNeighbours(a)) {
					if(aa.hasSpareValency()){
						continue atomLoop;
					}
				}
				a.setSpareValency(false);
			}
		}

		int svCount = 0;
		for(Atom a : atomCollection) {
			svCount += a.hasSpareValency() ? 1 :0;
		}

		/*
		 Reduce valency of atoms which cannot possibly have any of their bonds converted to double bonds
		 pick an atom which definitely does have spare valency to be the indicated hydrogen.
		*/
		Atom atomToReduceValencyAt =null;
		if (indicatedHydrogen!=null){
			if (getAtomByID(indicatedHydrogen)==null || !getAtomByID(indicatedHydrogen).hasSpareValency()){
				indicatedHydrogen = null;
			}
			else{
				atomToReduceValencyAt=getAtomByID(indicatedHydrogen);
			}
		}

		if((svCount % 2) == 1) {
			if (atomToReduceValencyAt == null){
				for(Atom a : atomCollection) {//try and find an atom with SV that neighbours only one atom with SV
					if(a.hasSpareValency()) {
						int atomsWithSV =0;
						for(Atom aa : getAtomNeighbours(a)) {
							if(aa.hasSpareValency()) {
								atomsWithSV++;
							}
						}
						if (atomsWithSV==1){
							atomToReduceValencyAt=a;
							break;
						}
					}
				}
				if (atomToReduceValencyAt==null){
					atomLoop: for(Atom a : atomCollection) {//try and find an atom with bridgehead atoms with SV on both sides c.f. phenoxastibinine ==10H-phenoxastibinine  else just pick the first atom with SV encountered
						if(a.hasSpareValency()) {
							if (atomToReduceValencyAt==null){
								atomToReduceValencyAt=a;
							}
							List<Atom> neighbours =getAtomNeighbours(a);
							if (neighbours.size()==2){
								for(Atom aa : neighbours) {
									if(getAtomNeighbours(aa).size() < 3){
										continue atomLoop;
									}
								}
								atomToReduceValencyAt=a;
								break;
							}
						}
					}
				}
			}
			atomToReduceValencyAt.setSpareValency(false);
			svCount--;
		}

		while(svCount > 0) {
			boolean foundTerminalFlag = false;
			boolean foundNonBridgeHeadFlag = false;
			boolean foundBridgeHeadFlag = false;
			for(Atom a : atomCollection) {
				if(a.hasSpareValency()) {
					int count = 0;
					for(Atom aa : getAtomNeighbours(a)) {
						if(aa.hasSpareValency()) {
							count++;
						}
					}
					if(count == 1) {
						for(Atom aa : getAtomNeighbours(a)) {
							if(aa.hasSpareValency()) {
								foundTerminalFlag = true;
								a.setSpareValency(false);
								aa.setSpareValency(false);
								findBondOrThrow(a, aa).addOrder(1);
								svCount -= 2;//Two atoms where for one of them this bond is the only double bond it can possible form
								break;
							}
						}
					}
				}
			}
			if(!foundTerminalFlag) {
				for(Atom a : atomCollection) {
					List<Atom> neighbours =getAtomNeighbours(a);
					if(a.hasSpareValency() && neighbours.size() < 3) {
						for(Atom aa : neighbours) {
							if(aa.hasSpareValency()) {
								foundNonBridgeHeadFlag = true;
								a.setSpareValency(false);
								aa.setSpareValency(false);
								findBondOrThrow(a, aa).addOrder(1);
								svCount -= 2;//Two atoms where one of them is not a bridge head
								break;
							}
						}
					}
					if(foundNonBridgeHeadFlag) break;
				}
				if(!foundNonBridgeHeadFlag){
					for(Atom a : atomCollection) {
						List<Atom> neighbours =getAtomNeighbours(a);
						if(a.hasSpareValency()) {
							for(Atom aa : neighbours) {
								if(aa.hasSpareValency()) {
									foundBridgeHeadFlag = true;
									a.setSpareValency(false);
									aa.setSpareValency(false);
									findBondOrThrow(a, aa).addOrder(1);
									svCount -= 2;//Two atoms where both of them are a bridge head e.g. necessary for something like coronene
									break;
								}
							}
						}
						if(foundBridgeHeadFlag) break;
					}
					if(!foundBridgeHeadFlag){
						throw new StructureBuildingException("Could not assign all higher order bonds.");
					}
				}
			}
		}
	}

	/**Gets a list of atoms in the fragment that connect to a specified atom
	 *
	 * @param atom The reference atom
	 * @return The list of atoms connected to the atom
	 * @throws StructureBuildingException
	 */
	List<Atom> getAtomNeighbours(Atom atom) throws StructureBuildingException {
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
	 * @return Incoming Valency
	 * @throws StructureBuildingException
	 */
	int getIntraFragmentIncomingValency(Atom atom) throws StructureBuildingException {
		int v = 0;
		for(Bond b :  atom.getBonds()) {
			//recalled atoms will be null if they are not part of this fragment
			if(b.getFromAtom() == atom) {
				Atom a =getAtomByID(b.getTo());
				if (a!=null && !a.getType().equals("suffix")){
					v += b.getOrder();
				}
			} else if(b.getToAtom() == atom) {
				Atom a =getAtomByID(b.getFrom());
				if (a!=null && !a.getType().equals("suffix")){
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
		for(Atom a : atomCollection) a.checkIncomingValency();
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

	int getDefaultInID() {
		return defaultInID;
	}

	void setDefaultInID(int defaultInID) {
		this.defaultInID=defaultInID;
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
	 * @param string
	 * @return
	 * @throws StructureBuildingException
	 */
	boolean hasLocant(String locant) throws StructureBuildingException {
		return getAtomByLocant(locant)!=null;
	}

	Integer getIndicatedHydrogen() {
		return indicatedHydrogen;
	}

	void setIndicatedHydrogen(Integer indicatedHydrogen) {
		this.indicatedHydrogen = indicatedHydrogen;
	}


	Atom getAtomByIdOrNextSuitableAtomOrThrow(int id, int additionalValencyRequired) throws StructureBuildingException {
		Atom a =getAtomByIdOrNextSuitableAtom(id, additionalValencyRequired, false);
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
	 * @param id
	 * @param additionalValencyRequired The increase in valency that will be required on the desired atom
	 * @return Atom
	 * @throws StructureBuildingException
	 */
	Atom getAtomByIdOrNextSuitableAtom(int id, int additionalValencyRequired, boolean takeIntoAccountOutValency) throws StructureBuildingException {
		List<Atom> atomList =getAtomList();
		int flag =0;
		Atom currentAtom=getAtomByIDOrThrow(id);
		int atomCounter=0;
		int atomListPosition=atomList.indexOf(currentAtom);
		int startingIndex =atomListPosition;

		do {//aromaticity preserved
			atomCounter++;
			if (atomListPosition >= atomList.size()){
				atomListPosition -=(atomList.size());
			}
			currentAtom=atomList.get(atomListPosition);
			if (atomCounter !=1 && currentAtom.getType().equals("suffix")){
				atomListPosition++;
				continue;
			}
			if (takeIntoAccountOutValency){
				if(ValencyChecker.checkValencyAvailableForBond(currentAtom, additionalValencyRequired + (currentAtom.hasSpareValency() ? 1 : 0) + currentAtom.getOutValency())){
					flag=1;
					break;
				}
			}
			else{
				if(ValencyChecker.checkValencyAvailableForBond(currentAtom, additionalValencyRequired + (currentAtom.hasSpareValency() ? 1 : 0))){
					flag=1;
					break;
				}
			}
			atomListPosition++;
		}
		while(atomCounter < atomList.size());

		atomListPosition =startingIndex;
		atomCounter=0;
		if (flag==0){

			do {//aromaticity dropped
				atomCounter++;
				if (atomListPosition >= atomList.size()){
					atomListPosition -=(atomList.size());
				}
				currentAtom=atomList.get(atomListPosition);
				if (atomCounter !=1 && currentAtom.getType().equals("suffix")){
					atomListPosition++;
					continue;
				}
				if (takeIntoAccountOutValency){
					if(ValencyChecker.checkValencyAvailableForBond(currentAtom, additionalValencyRequired + currentAtom.getOutValency())){
						flag=1;
						break;
					}
				}
				else{
					if(ValencyChecker.checkValencyAvailableForBond(currentAtom, additionalValencyRequired)){
						flag=1;
						break;
					}
				}
				atomListPosition++;
			}
			while(atomCounter < atomList.size());
		}
		if (flag==0){
			currentAtom=null;
		}
		return currentAtom;
	}

	/**
	 * Returns the id of the first atom that was added to this fragment
	 * @return
	 * @throws StructureBuildingException
	 */
	int getIdOfFirstAtom() throws StructureBuildingException {
		return getFirstAtom().getID();
	}

	/**
	 * Returns the the first atom that was added to this fragment
	 * @return
	 * @throws StructureBuildingException
	 */
	Atom getFirstAtom() throws StructureBuildingException {
		for (Atom a: atomCollection) {
			return a;
		}
		throw new StructureBuildingException("Fragment is empty");
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



