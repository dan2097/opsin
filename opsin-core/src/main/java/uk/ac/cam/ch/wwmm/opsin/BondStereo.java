package uk.ac.cam.ch.wwmm.opsin;

/**
 * Holds information about the positions of 2 atoms relative to a double bond allowing the specification of cis/trans stereochemistry
 * @author dl387
 *
 */
class BondStereo {

	private Atom[] atomRefs4;
	private BondStereoValue bondStereoValue;
	
	/**
	 * Possible values for a bondStereo element
	 * @author dl387
	 *
	 */
	enum BondStereoValue{
		CIS("C"),
		TRANS("T");

		private final String value;  
		BondStereoValue(String value){
			this.value = value;
		}
		@Override
		public String toString() {
			return value;
		}
	}

	/**
	 * Create a bondStereo from an array of 4 atoms. The 2nd and 3rd atoms of this array are connected via a double bond.
	 * The 1st and 4th atoms are at either end of this bond and indication is given as to whether they are cis or trans to each other.
	 * @param atomRefs4
	 * @param cOrT
	 */
	BondStereo(Atom[] atomRefs4, BondStereoValue cOrT) {
		if (atomRefs4.length !=4){
			throw new IllegalArgumentException("atomRefs4 must contain references to 4 atoms");
		}
		this.atomRefs4 = atomRefs4;
		this.bondStereoValue = cOrT;
	}
	
	Atom[] getAtomRefs4() {
		return atomRefs4;
	}
	void setAtomRefs4(Atom[] atomRefs4) {
		this.atomRefs4 = atomRefs4;
	}
	BondStereoValue getBondStereoValue() {
		return bondStereoValue;
	}
	void setBondStereoValue(BondStereoValue bondStereoValue) {
		this.bondStereoValue = bondStereoValue;
	}
}
