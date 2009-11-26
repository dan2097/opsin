package uk.ac.cam.ch.wwmm.opsin;

class OutID {
	/** Struct for an outID. As expected holds the id.
	 * Also holds the order of the bond that will be created when it is used (valency) e.g.  Eg. chloro 1, oxo 2
	 * setExplicitly says whether the outID absolutely definitely refers to that id or not.
	 * e.g. propyl is stored as prop-1-yl with this set to false while prop-2-yl has it set to true
	 *
	 * A reference to the fragment that this id is in is held
	 *
	 * Optionally a locant may be specified for what the outId should connect to if it is convenient to store such information. This is used in ester formation
	 */
	int id;
	int valency;
	boolean setExplicitly;
	Fragment frag;
	String locant = null;

	OutID(int id, int valency, Boolean setExplicitly, Fragment frag) {
		this.id = id;
		this.valency = valency;
		this.setExplicitly = setExplicitly;
		this.frag = frag;
	}

	OutID(int id, int valency, Boolean setExplicitly, Fragment frag, String locant) {
		this.id = id;
		this.valency = valency;
		this.setExplicitly = setExplicitly;
		this.frag = frag;
		this.locant = null;
	}
}
