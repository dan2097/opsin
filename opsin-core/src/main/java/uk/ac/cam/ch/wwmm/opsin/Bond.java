package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;
import nu.xom.Attribute;
import nu.xom.Element;

/**A bond, between two atoms.
 *
 * @author ptc24
 * @author dl387
 *
 */
class Bond {
	/** The Atom the bond comes from */
	private final Atom from;
	/** The Atom the bond goes to */
	private final Atom to;
	/** The bond order */
	private int order;

	static enum SMILES_BOND_DIRECTION{
		RSLASH,
		LSLASH
	}
	/** If this bond was built from SMILES can be set to either RSLASH or LSLASH. Subsequently read to add a bondStereoElement
	 * null by default*/
	private SMILES_BOND_DIRECTION smilesBondDirection = null;

	/**
	 * Holds the bondStereo object associated with this bond
	 * null by default
	 */
	private BondStereo bondStereo = null;

	/**
	 * If the bond is a fusion bond this will in the fused ring numberer be populated with the rings that it connects
	 */
	private final List<Ring> fusedRings = new ArrayList<Ring>(2);

	/** DO NOT CALL DIRECTLY EXCEPT FOR TESTING
	 * Creates a new Bond.
	 *
	 * @param from The Atom the bond comes from.
	 * @param to The Atom the bond goes to.
	 * @param order The bond order.
	 */
	Bond(Atom from, Atom to, int order) {
		if (from == to){
			throw new IllegalArgumentException("Bonds must be made between different atoms");
		}
		if (order < 1 || order > 3){
			throw new IllegalArgumentException("Bond order must be 1, 2 or 3");
		}
		if (from == null){
			throw new IllegalArgumentException("From atom was null!");
		}
		if (to == null){
			throw new IllegalArgumentException("To atom was null!");
		}
		this.from = from;
		this.to = to;
		this.order = order;
	}

	List<Ring> getFusedRings() {
		return fusedRings;
	}

	void  addFusedRing(Ring ring) {
		if (fusedRings.size()<2) {
			fusedRings.add(ring);
		}
	}


	/**Produces a nu.xom.Element corresponding to a CML bond tag.
	 * Has attributes of atomRefs2 and order.
	 *
	 * @return The CML element.
	 */
	Element toCMLBond() {
		Element elem = new Element("bond", XmlDeclarations.CML_NAMESPACE);
		elem.addAttribute(new Attribute("id", "a" + Integer.toString(from.getID())
				+ "_a" + Integer.toString(to.getID())));
		elem.addAttribute(new Attribute("atomRefs2", "a" + Integer.toString(from.getID())
				+ " a" + Integer.toString(to.getID())));
		if (order==1){
			elem.addAttribute(new Attribute("order", "S"));
		}
		else if (order==2){
			elem.addAttribute(new Attribute("order", "D"));
		}
		else if (order==3){
			elem.addAttribute(new Attribute("order", "T"));
		}
		else{
			elem.addAttribute(new Attribute("order", "unknown"));
		}
		if (bondStereo!=null){
			elem.appendChild(bondStereo.toCML());
		}
		return elem;
	}

	/**
	 * Gets from ID
	 * @return ID
	 */
	int getFrom() {
		return from.getID();
	}

	/**
	 * Gets to ID
	 * @return ID
	 */
	int getTo() {
		return to.getID();
	}

	/**Gets order.
    * @return*/
	int getOrder() {
		return order;
	}

	/**Sets order.
    * @param order*/
	void setOrder(int order) {
		this.order = order;
	}

	/**
	 * Gets from Atom
	 * @return Atom
	 */
	Atom getFromAtom() {
		return from;
	}

	/**
	 * Gets to Atom
	 * @return Atom
	 */
	Atom getToAtom() {
		return to;
	}

	/**Adds to the bond order.
	 *
	 * @param o The value to be added to the bond order.
	 */
	void addOrder(int o) {
		order += o;
	}

	/**
	 * Returns either null or RSLASH or LSLASH
	 * @return
	 */
	SMILES_BOND_DIRECTION getSmilesStereochemistry() {
		return smilesBondDirection;
	}

	void setSmilesStereochemistry(SMILES_BOND_DIRECTION bondDirection) {
		this.smilesBondDirection = bondDirection;
	}

	BondStereo getBondStereo() {
		return bondStereo;
	}

	void setBondStereo(BondStereo bondStereo) {
		this.bondStereo = bondStereo;
	}

	void setBondStereoElement(Atom[] atomRefs4, BondStereoValue cOrT) {
		bondStereo = new BondStereo(atomRefs4, cOrT);
	}

	/**
	 * Returns the atom at the other end of the bond to given atom
	 * @param atom
	 * @return
	 */
	Atom getOtherAtom(Atom atom) {
		if (from ==atom){
			return to;
		}
		else if (to ==atom){
			return from;
		}
		else{
			return null;
		}
	}

	@Override
	public int hashCode() {
		return from.hashCode() + to.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		Bond other = (Bond) obj;
		
		if (from == other.from && 
				to == other.to){
			return true;
		}
		if (from == other.to && 
				to == other.from){
			return true;
		}
		
		return false;
	}
}
