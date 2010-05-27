package uk.ac.cam.ch.wwmm.opsin;

import nu.xom.Attribute;
import nu.xom.Element;

class AtomParity {
	/**
	 * A dummy hydrogen atom. Used to represent an implicit hydrogen that is attached to a tetrahedral stereocentre
	 */
	static final Atom hydrogen = new Atom("H");
	private Atom[] atomRefs4;
	private int parity;
	
	/**
	 * Create an atomParity from an array of 4 atoms and the parity of the chiral determinant
	 * @param atomRefs4
	 * @param parity
	 * @throws StructureBuildingException
	 */
	AtomParity(Atom[] atomRefs4, int parity) throws StructureBuildingException {
		if (atomRefs4.length !=4){
			throw new StructureBuildingException("atomRefs4 must contain references to 4 atoms");
		}
		this.atomRefs4 = atomRefs4;
		this.parity = parity;
	}

	/**
	 * Serialises this object to CML
	 * @return
	 */
	Element toCML() {
		Element atomParityElement = new Element(XmlDeclarations.ATOMPARITY_EL);
		String atomRefsString ="";
		for (Atom atom : atomRefs4) {
			if (atomRefsString.equals("")){
				atomRefsString += "a" + atom.getID();
			}
			else{
				atomRefsString += " a" + atom.getID();
			}
		}
		atomParityElement.addAttribute(new Attribute(XmlDeclarations.ATOMREFS4_ATR, atomRefsString));
		atomParityElement.appendChild(Integer.toString(parity));
		return atomParityElement;
	}
	
	Atom[] getAtomRefs4() {
		return atomRefs4;
	}
	void setAtomRefs4(Atom[] atomRefs4) {
		this.atomRefs4 = atomRefs4;
	}
	int getParity() {
		return parity;
	}
	void setParity(int parity) {
		this.parity = parity;
	}
}
