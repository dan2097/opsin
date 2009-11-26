package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;

/**
 * An atom. Carries information about which fragment it is in, and an ID
 * number and a list of bonds that it is involved. It may also have other information such as
 * whether it has "spare valencies" due to unsaturation, its' charge, locant labels, stereochemistry and notes
 *
 * @author ptc24/dl387
 *
 */
class Atom {

	/**The (unique over the molecule) ID of the atom.*/
	private int ID;

	/**The atomic symbol of the atom. */
	private String element;

	/**The locants that pertain to the atom.*/
	private final List<String> locants = new ArrayList<String>();

	/**The formal charge on the atom.*/
	private int charge = 0;

	/**
	 * Holds the atomParity tag to be associated with this atom upon output
	 * null by default
	 */
	private Element atomParityElement = null;

	/**The bonds that involve the atom*/
	private final List<Bond> bonds = new ArrayList<Bond>();

	/**A hashmap in which notes can be associated with the atom. */
	private HashMap<String,String> notes = new HashMap<String, String>();

	/**The fragment to which the atom belongs.*/
	private Fragment frag;

	/**The number of not-to-implicit-hydrogen bonds going to or from the atom,
	 * that are not currently assigned. For example, benzene could be temporarily represented
	 * by six singly-bonded atoms, each with a spare valency of one, and later converted into
	 * a fully-specified valence structure.*/
	private int spareValency = 0;

	/**The total bond order of all bonds that are expected to be used for inter fragment bonding
	 * e.g. in butan-2-ylidene this would be 2 for the atom at position 2 and 0 for the other 3 */
	private int outValency = 0;

	/** Hydrogen count in CML, used subsequently to set the atom's valency when it's neighbours are known
	 * null if hydrogens are implicit*/
	private Integer explicitHydrogens;

	/** Null by default or set by the lambda convention/deduced from CML hydrogenCount.*/
	private Integer valency;

	/**
	 * Takes same values as type in Fragment. Useful for discriminating suffix atoms from other atoms when a suffix is incorporated into another fragments
	 */
	private String type;

	/**
	 * Is this atom in a ring. Default false. Set by the CycleDetector.
	 * Double bonds are only converted to spareValency if atom is in a ring
	 * Some suffixes have different meanings if an atom is part of a ring or not c.g. cyclohexanal vs ethanal
	 */
	private boolean atomIsInACycle =false;

	private static final Pattern matchElementSymbolLocant =Pattern.compile("[A-Z][a-z]?'*");
	private static final Pattern matchAminoAcidStyleLocant =Pattern.compile("([A-Z][a-z]?)('*)(\\d+[a-z]?'*)");

	/** Builds an Atom from scratch.
	 * @param ID The ID number, unique to the atom in the molecule being built
	 * @param locants A list of locants, typically numerical but can also include something like N
	 * @param element The atomic symbol of the chemical element
	 * @param frag the Fragment to contain the Atom
	 * @throws StructureBuildingException
	 */
	Atom(int ID, List<String> locants, String element, Fragment frag) throws StructureBuildingException {
		this(ID, element, frag);
		for (String locant : locants) {
			addLocant(locant);
		}
	}

	/** Builds an Atom from scratch.
	 * @param ID The ID number, unique to the atom in the molecule being built
	 * @param locant A locant, typically numerical but could also be something like N
	 * @param element The atomic symbol of the chemical element
	 * @param frag the Fragment to contain the Atom
	 * @throws StructureBuildingException
	 */
	Atom(int ID, String locant, String element, Fragment frag) throws StructureBuildingException {
		this(ID, element, frag);
		addLocant(locant);
	}

	/** Builds an Atom from scratch.
	 * @param ID The ID number, unique to the atom in the molecule being built
	 * @param element The atomic symbol of the chemical element
	 * @param frag the Fragment to contain the Atom
	 * @throws StructureBuildingException
	 */
	Atom(int ID, String element, Fragment frag) throws StructureBuildingException {
		if (frag==null){
			throw new StructureBuildingException("Atom is not in a fragment!");
		}
		if (element==null){
			throw new StructureBuildingException("Atom does not have an element!");
		}
		this.frag = frag;
		this.ID = ID;
		this.element = element;
		this.type =frag.getType();
	}

	/**Builds an Atom from a CML Atom tag. Looks at elementType and formalCharge
	 * attributes, hydrogenCount, and label tags contained within. id attributes are ignored.
	 *
	 * @param ID The ID number, unique to the atom in the molecule being built
	 * @param cmlAtom The nu.xom.Element for the Atom in CML
	 * @param frag the Fragment to contain the Atom
	 */
	Atom(int ID, Element cmlAtom, Fragment frag) {
		this.frag = frag;
		this.ID = ID;
		element = cmlAtom.getAttributeValue("elementType");
		this.type =frag.getType();
		Elements cmlLocants = cmlAtom.getChildElements("label");
		for(int i=0;i<cmlLocants.size();i++)
			addLocant(cmlLocants.get(0).getAttributeValue("value"));
		String chargeStr = cmlAtom.getAttributeValue("formalCharge");
		if(chargeStr != null){
			charge = Integer.parseInt(chargeStr);
		}
		String hcStr = cmlAtom.getAttributeValue("hydrogenCount");
		if(hcStr != null) {
			explicitHydrogens = Integer.parseInt(hcStr);
		}
	}

	/**Produces a nu.xom.Element for a CML Atom tag, containing
	 * values for id, elementType and (if appropriate) formalCharge, atomParity and
	 * embedded label tags.
	 *
	 * @return nu.xom.Element for a CML Atom tag
	 */
	Element toCMLAtom() {
		Element elem = new Element("atom");
		elem.addAttribute(new Attribute("id", "a" + Integer.toString(ID)));
		elem.addAttribute(new Attribute("elementType", element));
		if(charge != 0)
			elem.addAttribute(new Attribute("formalCharge", Integer.toString(charge)));
		if (explicitHydrogens!=null && explicitHydrogens==0){//prevent adding of implicit hydrogen
			elem.addAttribute(new Attribute("hydrogenCount", "0"));
		}
		if(atomParityElement != null){
			elem.appendChild(atomParityElement);
		}
		for(String l : locants) {
			Element locant = new Element("label");
			locant.addAttribute(new Attribute("value", l));
			locant.addAttribute(new Attribute("dictRef", "cmlDict:locant" ));
			elem.appendChild(locant);
		}
		return elem;
	}

	/**
	 * Uses the default, maximum and current incoming valency of the molecule to determine
	 * the likely current valency of the atom.
	 * The structureBuilder will of at this stage already rejected structures that violate the maxValency of an atom
	 * @return
	 */
	Integer determineValency() {
		if (valency != null){
			return valency;
		}
		Integer defaultValency =ValencyChecker.getDefaultValency(element, charge);
		Integer[] possibleValencies =ValencyChecker.getPossibleValencies(element, charge);
		int currentValency =getIncomingValency();

		if (defaultValency !=null && currentValency <= defaultValency){
			return defaultValency;
		}
		if (possibleValencies!=null) {
			for (Integer possibleValency : possibleValencies) {
				if (currentValency <= possibleValency){
					return possibleValency;
				}
			}
		}
		return null;
	}

	/**Adds a locant to the Atom. Other locants are preserved.
	 * Also associates the locant with the atom in the parent fragments hash
	 *
	 * @param locant The new locant
	 */
	void addLocant(String locant) {
		locants.add(locant);
		frag.addMappingToAtomLocantMap(locant, this);
	}

	/**Replaces all existing locants with a new one.
	 *
	 * @param locant The new locant
	 */
	void replaceLocant(String locant) {
		clearLocants();
		addLocant(locant);
	}

	void removeLocant(String locantToRemove) {
		int locantArraySize =locants.size();
		for (int i = locantArraySize -1; i >=0 ; i--) {
			if (locants.get(i).equals(locantToRemove)){
				locants.remove(i);
				frag.removeMappingFromAtomLocantMap(locantToRemove);
			}
		}
	}

	/**Removes all locants from the Atom.
	 *
	 */
	void clearLocants() {
		for (String l : locants) {
			frag.removeMappingFromAtomLocantMap(l);
		}
		locants.clear();
	}

	/**
	 * Removes only elementSymbolLocants: e.g. N, S', Se
	 */
	void removeElementSymbolLocants() {
		for (int i = locants.size()-1; i >=0; i--) {
			String l =locants.get(i);
			if (matchElementSymbolLocant.matcher(l).matches()){
				frag.removeMappingFromAtomLocantMap(l);
				locants.remove(i);
			}
		}
	}

	/**Checks if the Atom has a given locant.
	 *
	 * @param locant The locant to test for
	 * @return true if it has, false if not
	 * @throws StructureBuildingException
	 */
	boolean hasLocant(String locant) throws StructureBuildingException {
		for(String l : locants) {
			if(l.equals(locant))
				return true;
		}
		Matcher m = matchAminoAcidStyleLocant.matcher(locant);
		if (m.matches()){//e.g. N'5
			if (element.equals(m.group(1))){//element symbol
				if (!m.group(2).equals("") && (!hasLocant(m.group(1) +m.group(2)))){//has primes
					return false;//must have exact locant e.g. N'
				}
				if (OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(this, m.group(3))!=null){
					return true;
				}
			}
		}
		return false;
	}

	/**Gets the first locant for the Atom. This may be the locant that was initially
	 * specified, or the most recent locant specified using replaceLocant, or first
	 * locant to be added since the last invocation of clearLocants.
	 *
	 * @return The locant, or null if there is no locant
	 */
	String getFirstLocant() {
		if(locants.size() == 0)
			return null;
		return locants.get(0);
	}

	/**Returns the array of locants containing all locants associated with the atom
	 *
	 * @return The list of locants (may be empty)
	 */
	List<String> getLocants() {
		return locants;
	}

	/**Returns the subset of the locants which are element symbol locants e.g. N, S', Se
	 *
	 * @return The list of locants (may be empty)
	 */
	List<String> getElementSymbolLocants() {
		ArrayList<String> elementSymbolLocants =new ArrayList<String>();
        for (String l : locants) {
            if (matchElementSymbolLocant.matcher(l).matches()) {
                elementSymbolLocants.add(l);
            }
        }
		return elementSymbolLocants;
	}

	void setFrag(Fragment f) {
		frag = f;
	}

	Fragment getFrag() {
		return frag;
	}

	/**Gets the ID of the atom.
	 *
	 * @return The ID of the atom
	 */
	int getID() {
		return ID;
	}

	/**Gets the atomic symbol corresponding to the element of the atom.
	 *
	 * @return The atomic symbol corresponding to the element of the atom
	 */
	String getElement() {
		return element;
	}

	/**Sets the atomic symbol corresponding to the element of the atom.
	 *
	 * @param elem The atomic symbol corresponding to the element of the atom
	 */
	void setElement(String elem) {
		element = elem;
	}

	/**Gets the formal charge on the atom.
	 *
	 * @return The formal charge on the atom
	 */
	int getCharge() {
		return charge;
	}

	/**Sets the formal charge on the atom.
	 *
	 * @param c The formal charge on the atom
	 */
	void setCharge(int c) {
		charge = c;
		valency=null;//valency will now be non-standard
	}

	/**Adds a bond to the atom
	 *
	 * @param b The bond to be added
	 */
	void addBond(Bond b) {
		bonds.add(b);
	}

	/**Removes a bond to the atom
	 *
	 * @param b The bond to be removed
	 * @throws StructureBuildingException
	 */
	void removeBond(Bond b) throws StructureBuildingException {
		if(!bonds.remove(b)){
			throw new StructureBuildingException("Atom is not associated with this bond!");
		}
	}

	void checkIncomingValency() throws StructureBuildingException {
		if(!ValencyChecker.checkValency(this)) throw new StructureBuildingException();
	}

	/**Calculates the number of bonds connecting to the atom, excluding bonds to implicit
	 * hydrogens. Double bonds count as
	 * two bonds, etc. Eg ethene - both C's have an incoming valency of 2.
	 *
	 * @return Incoming Valency
	 */
	int getIncomingValency() {
		int v = 0;
		for(Bond b : bonds) {
			v += b.getOrder();
		}
		return v;
	}

	/**Gets the spare valency on the atom.
	 *
	 * @return The spare valency
	 */
	int getSpareValency() {
		return spareValency;
	}

	/**Sets the spare valency on the atom.
	 *
	 * @param sv The spare valency
	 */
	void setSpareValency(int sv) {
		spareValency = sv;
	}

	/**Adds to the spare valency on the atom.
	 *
	 * @param sv The spare valency to be added
	 */
	void addSpareValency(int sv) {
		spareValency += sv;
	}

	/**Subtracts from the spare valency on the atom.
	 *
	 * @param sv The spare valency to be subtracted
	 * @throws StructureBuildingException
	 */
	void subtractSpareValency(int sv) throws StructureBuildingException {
		if(spareValency < 1) throw new StructureBuildingException("Atom with locants " +
				locants + " has no spare valency to reduce.");
		spareValency -= sv;
	}

	/**Gets the total bond order of the bonds expected to be created from this atom for inter fragment bonding
	 *
	 * @return The outValency
	 */
	int getOutValency() {
		return outValency;
	}

	/**Adds to the total bond order of the bonds expected to be created from this atom for inter fragment bonding
	 *
	 * @param outV The outValency to be added
	 */
	void addOutValency(int outV) {
		outValency += outV;
	}


	/**
	 * Uses the explicitHydrogen attribute and incoming valency of the atom
	 * to calculate it's current valency
	 */
	void calculateValencyFromConnectivity() {
		if (explicitHydrogens != null){
			valency =explicitHydrogens + getIncomingValency();
		}
	}

	/**
	 * Allows the storage of a string of information associated with this atom in a hash
	 * @param key
	 * @param value
	 */
	void setNote(String key, String value) {
		notes.put(key,value);
	}

	/**
	 * Allows the retrieval of a string of information associated with this atom from a hash
	 * @param key
	 */
	String getNote(String key) {
		return notes.get(key);
	}

	List<Bond> getBonds() {
		return bonds;
	}

	/**Gets a list of atoms that connect to the atom
	 *
	 * @return The list of atoms connected to the atom
	 */
	List<Atom> getAtomNeighbours() throws StructureBuildingException {
		List<Atom> results = new ArrayList<Atom>();
		for(Bond b :  bonds) {
			if(b.getFromAtom() != this) {
				results.add(b.getFromAtom());
			} else if(b.getToAtom()!= this) {
				results.add(b.getToAtom());
			}
			else{
				throw new StructureBuildingException("Bond goes from an atom to the same atom!");
			}
		}
		return results;
	}


	Integer getExplicitHydrogens() {
		return explicitHydrogens;
	}

	void setExplicitHydrogens(Integer explicitHydrogens) {
		this.explicitHydrogens = explicitHydrogens;
	}

	Integer getValency() {
		return valency;
	}

	void setValency(Integer valency) {
		this.valency = valency;
	}

	String getType() {
		return type;
	}

	void setType(String type) {
		this.type = type;
	}

	boolean getAtomIsInACycle() {
		return atomIsInACycle;
	}

	/**
	 * Sets whether atom is in a cycle, true if it is
	 * @param atomIsInACycle
	 */
	void setAtomIsInACycle(boolean atomIsInACycle) {
		this.atomIsInACycle = atomIsInACycle;
	}

	Element getAtomParityElement() {
		return atomParityElement;
	}

	void setAtomParityElement(Element atomParityElement) {
		this.atomParityElement =atomParityElement;
	}

	void setAtomParityElement(String atomRefs4, int atomParity) {
		atomParityElement = new Element("atomParity");
		atomParityElement.addAttribute(new Attribute("atomRefs4",atomRefs4));
		atomParityElement.appendChild(Integer.toString(atomParity));
	}

	HashMap<String, String> getNotes() {
		return notes;
	}

	void setNotes(HashMap<String, String> notes) {
		this.notes = notes;
	}

	/**
	 * Checks if the valency of this atom allows it to have the amount of spare valency that the atom currently has
	 * May reduce the spare valency on the atom to be consistant with the valency of the atom
	 * Sets a note if a nitrogen with 3 bonds and a spare valency is encountered e.g. N-methylpyridine
	 * Does nothing if the atom has no spare valency
	 * @param takeIntoAccountExternalBonds
	 * @throws StructureBuildingException
	 */
	void ensureSVIsConsistantWithValency(boolean takeIntoAccountExternalBonds) throws StructureBuildingException {
		if (spareValency > 0 ){
			Integer maxValency;
			if (valency!=null){
				maxValency=valency;
			}
			else{
				if (charge==0 && ValencyChecker.getHWValency(element)!=null){
					maxValency =ValencyChecker.getHWValency(element);
				}
				else{
					maxValency= ValencyChecker.getMaximumValency(element, charge);
				}
			}
			if (maxValency !=null){
				int maxSpareValency;
				if (takeIntoAccountExternalBonds){
					maxSpareValency =maxValency-getIncomingValency();
				}
				else{
					maxSpareValency =maxValency-frag.getIntraFragmentIncomingValency(this);
				}
				if (maxSpareValency < spareValency){
					if (element.equals("N")&& spareValency==1 && charge==0){//special case where a charge has erroneously been omitted e.g. 1-methylpyridine
						setNote("Possibly Should Be Charged", "1");
					}
					if (maxSpareValency > 0 ){
						setSpareValency(maxSpareValency);
					}
					else{
						setSpareValency(0);
					}
				}
			}
		}
	}
}
