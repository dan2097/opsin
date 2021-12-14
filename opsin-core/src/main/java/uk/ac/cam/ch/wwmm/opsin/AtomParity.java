package uk.ac.cam.ch.wwmm.opsin;

/**
 * Hold information about 4 atoms and their chiral determinant allowing the description of tetrahedral stereochemistry
 * @author dl387
 *
 */
class AtomParity {
	/**
	 * A dummy hydrogen atom. Used to represent an implicit hydrogen that is attached to a tetrahedral stereocentre
	 */
	static final Atom hydrogen = new Atom(ChemEl.H);

	/**
	 * A dummy hydrogen atom. Used to represent the hydrogen that replaced a hydroxy at a tetrahedral stereocentre
	 */
	static final Atom deoxyHydrogen = new Atom(ChemEl.H);
	
	/**
	 * Typical case where of an absolute stereocentre in group number 1
	 */
	private static final StereoGroup ABSOLUTE_STEREOGROUP = new StereoGroup(StereoGroupType.Abs);
	

	private Atom[] atomRefs4;
	private int parity;
	private StereoGroup stereoGroup = ABSOLUTE_STEREOGROUP;
	
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

	void setStereoGroup(StereoGroup stroGrp) {
		this.stereoGroup = stroGrp;
	}

	StereoGroup getStereoGroup() {
		return this.stereoGroup;
	}
}
