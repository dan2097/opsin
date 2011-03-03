package uk.ac.cam.ch.wwmm.opsin;

/**
 * Struct for a FunctionalAtom. As expected holds the atom.
 * This is used to indicate, for example, that this atom may form an ester
 * 
 * @author dl387
 *
 */
class FunctionalAtom {
	private final Atom atom;

	FunctionalAtom(Atom atom) {
		this.atom = atom;
	}

	Atom getAtom() {
		return atom;
	}
}
