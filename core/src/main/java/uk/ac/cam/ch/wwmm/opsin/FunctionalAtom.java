package uk.ac.cam.ch.wwmm.opsin;

class FunctionalAtom {
	/** Struct for a FunctionalAtom. As expected holds the atom.
	 */
	private Atom atom;

	FunctionalAtom(Atom atom) {
		this.atom = atom;
	}

	Atom getAtom() {
		return atom;
	}

	void setAtom(Atom atom) {
		this.atom = atom;
	}
}
