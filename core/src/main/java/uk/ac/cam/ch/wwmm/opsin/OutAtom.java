package uk.ac.cam.ch.wwmm.opsin;

class OutAtom {
	/** Struct for an OutAtom. As expected holds a reference to an atom.
	 * However if setExplicitly is not true then it is not true then it is not absolutely definitely this amount that is referred to
	 * e.g. propyl is stored as prop-1-yl with this set to false while prop-2-yl has it set to true
	 * 
	 * Also holds the order of the bond that will be created when it is used (valency) e.g.  Eg. chloro 1, oxo 2
	 *
	 * Optionally a locant may be specified for what the outAtom should connect to if it is convenient to store such information. This is used in ester formation
	 */
	private Atom atom;
	private int valency;
	private boolean setExplicitly;
	private String locant;
	
	OutAtom(Atom atom, int valency, Boolean setExplicitly) {
		this(atom, valency, setExplicitly, null);
	}

	OutAtom(Atom atom, int valency, Boolean setExplicitly, String locant) {
		this.atom = atom;
		this.valency = valency;
		this.setExplicitly = setExplicitly;
		this.locant = null;
		if (setExplicitly){
			atom.addOutValency(valency);
		}
	}
	Atom getAtom() {
		return atom;
	}

	void setAtom(Atom atom) {
		this.atom = atom;
	}

	int getValency() {
		return valency;
	}

	void setValency(int valency) {
		if (setExplicitly){
			atom.addOutValency(valency -this.valency);
		}
		this.valency = valency;
	}

	boolean isSetExplicitly() {
		return setExplicitly;
	}

	void setSetExplicitly(boolean setExplicitly) {
		if (!this.setExplicitly && setExplicitly){
			atom.addOutValency(valency);
		}
		else if (this.setExplicitly && !setExplicitly){
			atom.addOutValency(-valency);
		}
		this.setExplicitly = setExplicitly;
	}

	String getLocant() {
		return locant;
	}

	void setLocant(String locant) {
		this.locant = locant;
	}
}
