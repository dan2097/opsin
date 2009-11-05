package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;

import nu.xom.Attribute;
import nu.xom.Element;

/**A bond, between two atoms.
 *
 * @author ptc24/dl387
 *
 */
class Bond {
	/** The Atom the bond comes from */
	private Atom from;
	/** The Atom the bond goes to */
	private Atom to;
	/** The bond order */
	private int order;

	/**If the bond can be either cis or trans this can be defined by setting this to either "E" or "Z"
	 * "e" or "z" indicate relative stereochemistry
	 * null by default*/
	private String stereochemistry = null;
	
	/**
	 * Holds the bondStereo tag to be associated with this bond upon output
	 * null by default
	 */
	private Element bondStereoElement = null;
	
	private int smilesNumber;
	
	/**
	 * it the bond is fusion bond it contains the rings that it connects
	 */
	private ArrayList<Ring> fusedRings = new ArrayList<Ring>(2);

	/**Creates a new Bond.
	 *
	 * @param from The Atom the bond comes from.
	 * @param to The Atom the bond goes to.
	 * @param order The bond order.
	 */
	Bond(Atom from, Atom to, int order) {
		this.from = from;
		this.to = to;
		this.order = order;
		smilesNumber = 0;
	}
	
	public ArrayList<Ring> getFusedRings() {
		return fusedRings;
	}
	
	public void  addFusedRing(Ring ring) {
		if (fusedRings.size()<2) fusedRings.add(ring);
	}
	

	/**Produces a nu.xom.Element corresponding to a CML bond tag.
	 * Has attributes of atomRefs2 and order.
	 *
	 * @return The CML element.
	 */
	Element toCMLBond() {
		Element elem = new Element("bond");
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
		if (bondStereoElement!=null){
			bondStereoElement.addAttribute(new Attribute("convention","cmlDict:cistrans"));
			elem.appendChild(bondStereoElement);
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

	/**Gets order.*/
	int getOrder() {
		return order;
	}

	/**Sets order.*/
	void setOrder(int o) {
		order = o;
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
	 * Returns either null or one of E,Z,e,z
	 * @return
	 */
	String getStereochemistry() {
		return stereochemistry;
	}

	void setStereochemistry(String stereochemistry) {
		this.stereochemistry = stereochemistry;
	}

	Element getBondStereoElement() {
		return bondStereoElement;
	}

	void setBondStereoElement(Element bondStereoElement) {
		this.bondStereoElement = bondStereoElement;
	}

	void setBondStereoElement(String atomRefs4, String CorT) {
		bondStereoElement = new Element("bondStereo");
		bondStereoElement.addAttribute(new Attribute("atomRefs4", atomRefs4));
		bondStereoElement.appendChild(CorT);
	}

	int getSMILESNumber() {
		return smilesNumber;
	}

	int assignSMILESNumber(IDManager idm) {
		smilesNumber = idm.getNextID();
		return smilesNumber;
	}


}
