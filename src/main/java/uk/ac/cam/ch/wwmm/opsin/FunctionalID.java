package uk.ac.cam.ch.wwmm.opsin;

public class FunctionalID {
	/** Struct for a FunctionalID. As expected holds the id.
	 * A reference to the fragment that this id is in is held
	 */
	int id;
	Fragment frag;

	FunctionalID(int id, Fragment frag) {
		this.id = id;
		this.frag = frag;
	}
}
