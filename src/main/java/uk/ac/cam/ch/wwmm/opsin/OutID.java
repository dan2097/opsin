package uk.ac.cam.ch.wwmm.opsin;



public class OutID {
	/** Struct for an outID. As expected holds the id.
	 * Also holds the order of the bond that will be created when it is used (valency) e.g.  Eg. chloro 1, oxo 2
	 * setExplicitly says whether the outID absolutely definitely refers to that id or not.
	 * e.g. propyl is stored as prop-1-yl with this set to false while prop-2-yl has it set to true
	 * 
	 * A reference to the fragment that this id is in is held
	 * Optionally a locant that the outID is expected to go to is stored. This will be used for names like 3,4-methylenedioxyphen-1-yl (the outIDs on the oxys would have this set)
	 */
	int id;
	int valency;
	boolean setExplicitly;
	Fragment frag;
	String locant;
	BuildResults buildResults = null;//TODO get rid of this abomination
	
	void setBuildResults(BuildResults buildResults) {
		this.buildResults = buildResults;
	}

	OutID(int id, int valency, Boolean setExplicitly, Fragment frag, String locant) {
		this.id = id;
		this.valency = valency;
		this.setExplicitly = setExplicitly;
		this.frag = frag;
		this.locant = locant;
	}
}
