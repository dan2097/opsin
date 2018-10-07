package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**
 * An atom. Carries information about which fragment it is in, and an ID
 * number and a list of bonds that it is involved. It may also have other information such as
 * whether it has "spare valencies" due to unsaturation, its charge, locant labels, stereochemistry and notes
 *
 * @author ptc24
 * @author dl387
 *
 */
class Atom {

	/**The (unique over the molecule) ID of the atom.*/
	private final int id;

	/**The chemical element of the atom. */
	private ChemEl chemEl;

	/**The locants that pertain to the atom.*/
	private final List<String> locants = new ArrayList<String>(2);

	/**The formal charge on the atom.*/
	private int charge = 0;
	
	/**The isotope of the atom. Null if not defined explicitly.*/
	private Integer isotope = null;

	/**
	 * Holds the atomParity object associated with this object
	 * null by default
	 */
	private AtomParity atomParity = null;

	/**The bonds that involve the atom*/
	private final List<Bond> bonds = new ArrayList<Bond>(4);

	/**A map between PropertyKey s as declared here and useful atom properties, usually relating to some kind of special case. */
	@SuppressWarnings("rawtypes")
	private final Map<PropertyKey, Object> properties = new HashMap<PropertyKey, Object>();
	/** A set of atoms that were equally plausible to perform functional replacement on */
	static final PropertyKey<Set<Atom>> AMBIGUOUS_ELEMENT_ASSIGNMENT = new PropertyKey<Set<Atom>>("ambiguousElementAssignment");
	/** The atom class which will be output when serialised to SMILES. Useful for distinguishing attachment points */
	static final PropertyKey<Integer> ATOM_CLASS = new PropertyKey<Integer>("atomClass");
	/** Used on wildcard atoms to indicate their meaning */
	static final PropertyKey<String> HOMOLOGY_GROUP = new PropertyKey<String>("homologyGroup");
	/** Used on wildcard atoms to indicate that they are a position variation bond */
	static final PropertyKey<List<Atom>> POSITION_VARIATION_BOND = new PropertyKey<List<Atom>>("positionVariationBond");
	/** The hydrogen count as set in the SMILES*/
	static final PropertyKey<Integer> SMILES_HYDROGEN_COUNT = new PropertyKey<Integer>("smilesHydrogenCount");
	/** The oxidation number as specified by Roman numerals in the name*/
	static final PropertyKey<Integer> OXIDATION_NUMBER = new PropertyKey<Integer>("oxidationNumber");
	/** Is this atom the carbon of an aldehyde? (however NOT formaldehyde)*/
	static final PropertyKey<Boolean> ISALDEHYDE = new PropertyKey<Boolean>("isAldehyde");
	/** Indicates that this atom is an anomeric atom in a cyclised carbohydrate*/
	static final PropertyKey<Boolean> ISANOMERIC = new PropertyKey<Boolean>("isAnomeric");
	/** Transient integer used to indicate traversal of fragments*/
	static final PropertyKey<Integer> VISITED = new PropertyKey<Integer>("visited");

	/**The fragment to which the atom belongs.*/
	private Fragment frag;

	/** Whether an atom is part of a delocalised set of double bonds. A double bond in a kekule structure
	 * can be mapped to a single bond with this attribute set to true on both atoms that were in the double bond
	 * For example, benzene could be temporarily represented by six singly-bonded atoms, each with a set
	 * spare valency attribute , and later converted into a fully-specified valence structure.*/
	private boolean spareValency = false;

	/**The total bond order of all bonds that are expected to be used for inter fragment bonding
	 * e.g. in butan-2-ylidene this would be 2 for the atom at position 2 and 0 for the other 3 */
	private int outValency = 0;

	/** Null by default or set by the lambda convention.*/
	private Integer lambdaConventionValency;
	
	/** Null by default or set by the SMILES builder*/
	private Integer minimumValency;
	
	/** Can this atom have implicit hydrogen? True unless explicitly set otherwise otherwise*/
	private boolean implicitHydrogenAllowed = true;

	/** This is modified by ium/ide/ylium/uide and is used to choose the appropriate valency for the atom*/
	private int protonsExplicitlyAddedOrRemoved = 0;

	/**
	 * Takes same values as type in Fragment. Useful for discriminating suffix atoms from other atoms when a suffix is incorporated into another fragments
	 */
	private String type;

	/**
	 * Is this atom in a ring. Default false. Set by the CycleDetector.
	 * Double bonds are only converted to spareValency if atom is in a ring
	 * Some suffixes have different meanings if an atom is part of a ring or not c.g. cyclohexanal vs ethanal
	 */
	private boolean atomIsInACycle = false;

	/**
	 * Builds an Atom from scratch.
	 * GENERALLY EXCEPT FOR TESTING SHOULD NOT BE CALLED EXCEPT FROM THE FRAGMANAGER
	 * @param id The ID number, unique to the atom in the molecule being built
	 * @param chemlEl The chemical element
	 * @param frag the Fragment to contain the Atom
	 */
	Atom(int id, ChemEl chemlEl, Fragment frag) {
		if (frag == null){
			throw new IllegalArgumentException("Atom is not in a fragment!");
		}
		if (chemlEl == null){
			throw new IllegalArgumentException("Atom does not have an element!");
		}
		this.frag = frag;
		this.id = id;
		this.chemEl = chemlEl;
		this.type =frag.getType();
	}
	
	/** Used to build a DUMMY atom.
	 * Does not have an id/frag/type as would be expected for a proper atom
	 * @param chemlEl The chemical element
	 */
	Atom(ChemEl chemlEl){
		this.chemEl = chemlEl;
		this.id = 0;
	}
	
	/**
	 * Uses the lambdaConventionValency or if that is not available
	 * the default valency assuming this is >= the current valency
	 * If not then allowed the chemically sensible valencies of the atom are checked with the one that is closest and >= to the current valency
	 * being returned. If the valency has still not been determined the current valency i.e. assuming the atom to have 0 implicit hydrogen is returned.
	 * This is the correct behaviour for inorganics. For p block elements it means that OPSIN does not believe the atom to be in a valid valency (too high)
	 * 
	 * if considerOutValency is true, the valency that will be used to form bonds using the outAtoms is
	 * taken into account i.e. if any radicals were used to form bonds
	 * @param considerOutValency
	 * @return
	 */
	int determineValency(boolean considerOutValency) {
		if (lambdaConventionValency != null){
			return lambdaConventionValency + protonsExplicitlyAddedOrRemoved;
		}
		int currentValency = getIncomingValency();
		if (considerOutValency){
			currentValency += outValency;
		}
		Integer calculatedMinValency = minimumValency == null ? null : minimumValency + protonsExplicitlyAddedOrRemoved;
		if (charge ==0 || protonsExplicitlyAddedOrRemoved != 0){
			Integer defaultValency = ValencyChecker.getDefaultValency(chemEl);
			if (defaultValency != null){
				defaultValency += protonsExplicitlyAddedOrRemoved;
				if (currentValency <= defaultValency && (calculatedMinValency == null || defaultValency >= calculatedMinValency)){
					return defaultValency;
				}
			}
		}
		Integer[] possibleValencies = ValencyChecker.getPossibleValencies(chemEl, charge);
		if (possibleValencies != null) {
			if (calculatedMinValency != null  && calculatedMinValency >= currentValency){
				return calculatedMinValency;
			}
			for (Integer possibleValency : possibleValencies) {
				if (calculatedMinValency != null && possibleValency < calculatedMinValency){
					continue;
				}
				if (currentValency <= possibleValency){
					return possibleValency;
				}
			}
		}
		if (calculatedMinValency != null && calculatedMinValency >= currentValency){
			return calculatedMinValency;
		}
		else{
			return currentValency;
		}
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
	void replaceLocants(String locant) {
		clearLocants();
		addLocant(locant);
	}

	void removeLocant(String locantToRemove) {
		int locantArraySize = locants.size();
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
		for (int i = 0, l = locants.size(); i < l; i++) {
			frag.removeMappingFromAtomLocantMap(locants.get(i));
		}
		locants.clear();
	}

	/**
	 * Removes only elementSymbolLocants: e.g. N, S', Se
	 */
	void removeElementSymbolLocants() {
		for (int i = locants.size() - 1; i >= 0; i--) {
			String locant = locants.get(i);
			if (MATCH_ELEMENT_SYMBOL_LOCANT.matcher(locant).matches()){
				frag.removeMappingFromAtomLocantMap(locant);
				locants.remove(i);
			}
		}
	}
	
	/**
	 * Removes all locants other than elementSymbolLocants (e.g. N, S', Se)
	 * Hence removes numeric locants and greek locants
	 */
	void removeLocantsOtherThanElementSymbolLocants() {
		for (int i = locants.size() - 1; i >= 0; i--) {
			String locant = locants.get(i);
			if (!MATCH_ELEMENT_SYMBOL_LOCANT.matcher(locant).matches()){
				frag.removeMappingFromAtomLocantMap(locant);
				locants.remove(i);
			}
		}
	}

	/**Checks if the Atom has a given locant.
	 *
	 * @param locant The locant to test for
	 * @return true if it has, false if not
	 */
	boolean hasLocant(String locant) {
		if (locants.contains(locant)) {
			return true;
		}
		Matcher m = MATCH_AMINOACID_STYLE_LOCANT.matcher(locant);
		if (m.matches()){//e.g. N'5
			if (chemEl.toString().equals(m.group(1))){//element symbol
				if (!m.group(2).equals("") && (!hasLocant(m.group(1) +m.group(2)))){//has primes
					return false;//must have exact locant e.g. N'
				}
				if (OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(this, m.group(3)) != null){
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
		return locants.size() > 0 ? locants.get(0) : null;
	}

	/**Returns the array of locants containing all locants associated with the atom
	 *
	 * @return The list of locants (may be empty)
	 */
	List<String> getLocants() {
		return Collections.unmodifiableList(locants);
	}

	/**Returns the subset of the locants which are element symbol locants e.g. N, S', Se
	 *
	 * @return The list of locants (may be empty)
	 */
	List<String> getElementSymbolLocants() {
		List<String> elementSymbolLocants = new ArrayList<String>(1);
		for (int i = 0, l = locants.size(); i < l; i++) {
			String locant = locants.get(i);
            if (MATCH_ELEMENT_SYMBOL_LOCANT.matcher(locant).matches()) {
                elementSymbolLocants.add(locant);
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
		return id;
	}

	/**Gets the chemical element corresponding to the element of the atom.
	 *
	 * @return The chemical element corresponding to the element of the atom
	 */
	ChemEl getElement() {
		return chemEl;
	}

	/**Sets the chemical element corresponding to the element of the atom.
	 *
	 * @param chemEl The chemical element corresponding to the element of the atom
	 */
	void setElement(ChemEl chemEl) {
		this.chemEl = chemEl;
	}

	/**Gets the formal charge on the atom.
	 *
	 * @return The formal charge on the atom
	 */
	int getCharge() {
		return charge;
	}
	
	/**Modifies the charge of this atom by the amount given. This can be any integer
	 * The number of protons changed is noted so as to calculate the correct valency for the atom. This can be any integer.
	 * For example ide is the loss of a proton so is charge=-1, protons =-1
	 * @param charge
	 * @param protons
	 */
	void addChargeAndProtons(int charge, int protons){
		this.charge += charge;
		protonsExplicitlyAddedOrRemoved+=protons;
	}

	/**Sets the formal charge on the atom.
	 * NOTE: make sure to update protonsExplicitlyAddedOrRemoved if necessary
	 * 
	 * @param c The formal charge on the atom
	 */
	void setCharge(int c) {
		charge = c;
	}
	
	 /** 
	 * Sets the formal charge and number of protonsExplicitlyAddedOrRemoved to 0
	 */
	void neutraliseCharge() {
		charge = 0;
		protonsExplicitlyAddedOrRemoved = 0;
	}

	/**
	 * Gets the mass number of the atom or null if not explicitly defined
	 * e.g. 3 for tritium
	 * @return
	 */
	Integer getIsotope() {
		return isotope;
	}

	/**
	 * Sets the mass number of the atom explicitly
	 * @param isotope
	 */
	void setIsotope(Integer isotope) {
		if (isotope != null && isotope < chemEl.ATOMIC_NUM) {
			throw new RuntimeException("Isotopic mass cannot be less than the element's number of protons: " + chemEl.toString() + " " + isotope + " < " + chemEl.ATOMIC_NUM );
		}
		this.isotope = isotope;
	}

	/**Adds a bond to the atom
	 *
	 * @param b The bond to be added
	 */
	void addBond(Bond b) {
		if (bonds.contains(b)){
			throw new IllegalArgumentException("Atom already has given bond (This is not allowed as this would give two bonds between the same atoms!)");
		}
		bonds.add(b);
	}

	/**Removes a bond to the atom
	 *
	 * @param b The bond to be removed
     * @return whether bond was present
	 */
	boolean removeBond(Bond b) {
		return bonds.remove(b);
	}

	/**Calculates the number of bonds connecting to the atom, excluding bonds to implicit
	 * hydrogens. Double bonds count as
	 * two bonds, etc. Eg ethene - both C's have an incoming valency of 2.
	 *
	 * @return Incoming Valency
	 */
	int getIncomingValency() {
		int v = 0;
		for (int i = 0, l = bonds.size(); i < l; i++) {
			v += bonds.get(i).getOrder();
		}
		return v;
	}

	int getProtonsExplicitlyAddedOrRemoved() {
		return protonsExplicitlyAddedOrRemoved;
	}

	void setProtonsExplicitlyAddedOrRemoved(int protonsExplicitlyAddedOrRemoved) {
		this.protonsExplicitlyAddedOrRemoved = protonsExplicitlyAddedOrRemoved;
	}
	
	/**Does the atom have spare valency to form double bonds?
	 *
	 * @return true if atom has spare valency
	 */
	boolean hasSpareValency() {
		return spareValency;
	}

	/**Set whether an atom has spare valency
	 *
	 * @param sv The spare valency
	 */
	void setSpareValency(boolean sv) {
		spareValency = sv;
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

	List<Bond> getBonds() {
		return Collections.unmodifiableList(bonds);
	}
	
	int getBondCount() {
		return bonds.size();
	}

	/**Gets a list of atoms that connect to the atom
	 *
	 * @return The list of atoms connected to the atom
	 */
	List<Atom> getAtomNeighbours(){
		int bondCount = bonds.size();
		List<Atom> results = new ArrayList<Atom>(bondCount);
		for (int i = 0; i < bondCount; i++) {
			results.add(bonds.get(i).getOtherAtom(this));
		}
		return results;
	}

	Integer getLambdaConventionValency() {
		return lambdaConventionValency;
	}

	void setLambdaConventionValency(Integer valency) {
		this.lambdaConventionValency = valency;
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

	AtomParity getAtomParity() {
		return atomParity;
	}

	void setAtomParity(AtomParity atomParity) {
		this.atomParity = atomParity;
	}

	void setAtomParity(Atom[] atomRefs4, int parity) {
		atomParity = new AtomParity(atomRefs4, parity);
	}
	
	Integer getMinimumValency() {
		return minimumValency;
	}

	void setMinimumValency(Integer minimumValency) {
		this.minimumValency = minimumValency;
	}
	
	boolean getImplicitHydrogenAllowed() {
		return implicitHydrogenAllowed;
	}

	void setImplicitHydrogenAllowed(boolean implicitHydrogenAllowed) {
		this.implicitHydrogenAllowed = implicitHydrogenAllowed;
	}

    @SuppressWarnings("unchecked")
	<T> T getProperty(PropertyKey<T> propertyKey) {
        return (T) properties.get(propertyKey);
    }

	<T> void setProperty(PropertyKey<T> propertyKey, T value) {
		properties.put(propertyKey, value);
	}

	/**
	 * Checks if the valency of this atom allows it to have the amount of spare valency that the atom currently has
	 * May reduce the spare valency on the atom to be consistent with the valency of the atom
	 * Does nothing if the atom has no spare valency
	 * @param takeIntoAccountExternalBonds
	 * @throws StructureBuildingException
	 */
	void ensureSVIsConsistantWithValency(boolean takeIntoAccountExternalBonds) throws StructureBuildingException {
		if (spareValency) {
			Integer maxValency;
			if (lambdaConventionValency != null) {
				maxValency = lambdaConventionValency + protonsExplicitlyAddedOrRemoved;
			}
			else{
				Integer hwValency = ValencyChecker.getHWValency(chemEl);
				if (hwValency == null) {
					throw new StructureBuildingException(chemEl + " is not expected to be aromatic!");
				}
				maxValency = hwValency + protonsExplicitlyAddedOrRemoved;
			}
			int maxSpareValency;
			if (takeIntoAccountExternalBonds) {
				maxSpareValency = maxValency - getIncomingValency() - outValency;
			}
			else{
				maxSpareValency = maxValency - frag.getIntraFragmentIncomingValency(this);
			}
			if (maxSpareValency < 1) {
				setSpareValency(false);
			}
		}
	}
	
	/**
	 * Returns the the first bond in the atom's bond list or null if it has no bonds
	 * @return
	 */
	Bond getFirstBond() {
		if (bonds.size() > 0){
			return bonds.get(0);
		}
		return null;
	}
	
	/**Gets the bond between this atom and a given atom
	 *
	 * @param a The atom to find a bond to
	 * @return The bond, or null if there is no bond
	 */
	Bond getBondToAtom(Atom a) {
		for (int i = 0, l = bonds.size(); i < l; i++) {
			Bond b = bonds.get(i);
			if(b.getOtherAtom(this) == a){
				return b;
			}
		}
		return null;
	}

	/**Gets the bond between this atom and a given atom, throwing if fails.
	 *
	 * @param a The atom to find a bond to
	 * @return The bond found
	 * @throws StructureBuildingException
	 */
	Bond getBondToAtomOrThrow(Atom a) throws StructureBuildingException {
		Bond b = getBondToAtom(a);
		if(b == null){
			throw new StructureBuildingException("Couldn't find specified bond");
		}
		return b;
	}
}
