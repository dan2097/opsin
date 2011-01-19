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
	 */
	AtomParity(Atom[] atomRefs4, int parity){
		if (atomRefs4.length !=4){
			throw new IllegalArgumentException("atomRefs4 must contain references to 4 atoms");
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
		StringBuilder atomRefsSb = new StringBuilder();
		for(int i=0; i<atomRefs4.length-1; i++) {
			atomRefsSb.append('a');
			atomRefsSb.append(atomRefs4[i].getID());
			atomRefsSb.append(' ');
		}
		atomRefsSb.append('a');
		atomRefsSb.append(atomRefs4[atomRefs4.length-1].getID());
		atomParityElement.addAttribute(new Attribute(XmlDeclarations.ATOMREFS4_ATR, atomRefsSb.toString()));
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
