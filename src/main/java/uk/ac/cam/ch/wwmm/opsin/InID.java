package uk.ac.cam.ch.wwmm.opsin;

public class InID {
	/** Struct for an InID. As expected holds the id.
	 * valency gives the order of the bond to be created when connected. Important for suffixes
	 * A reference to the fragment that this id is in is held
	 */
	int id;
	int valency;
	Fragment frag;

	InID(int id, int valency, Fragment frag) {
		this.id = id;
		this.valency = valency;
		this.frag = frag;
	}
}